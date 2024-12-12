function [MeshGeometry] = axisymmetricGeometry(Specimen,Indenter)
%% Specimen parameters
    R = Specimen.R; % Specimen radius
    H = Specimen.H; % Specimen height
    NR = Specimen.NR; % Number of partitions along radial direction
    NH = Specimen.NH; % Number of partitions along axial direction
    
    R_bias = Specimen.R_bias; % Radial partition bias - X=x*|x/L|^R_bias
    H_bias = Specimen.H_bias; % Axial partition bias - Z=z*|z/H|^H_bias
    benchmark_flag = Specimen.benchmark_flag; % benchmark flag
    elementType = Specimen.elementType; % hex20 / hex8
    ignore_formula = Specimen.ignore_formula; % flag for calculating alpha and rho
    
    % Calculate the angle of the cylindrical sector, alpha.
    if ignore_formula % use default nodal position scheme (see eqns. (1.1-1.2) and Fig. 1b in paper)
        alpha = Specimen.alpha;
        % make sure that alpha is a divisor of 360 so that the sector could be 
        % replicated N_repeat times if a full 3D model is needed, such as for benchmarking
        % purposes. 
        N_repeat_trial = 360/alpha; % number of repititions of sectors needed to complete an entire 3D cylinder
        if N_repeat_trial==floor(N_repeat_trial) % 
            N_repeat = N_repeat_trial;
        else % Update alpha if needed
            N_repeat = round(N_repeat_trial);
            alpha = 360/N_repeat; %actual angle
            warning('Alpha angle was changed to %d',alpha); % notify the user if alpha was changed
        end
    else
        rho = 0.4; % radius of pentahedral elements 
        alpha = rad2deg(2*(atan((R/rho)^(1/(NR-1)))-pi/4));% trial angle
        N_repeat_trial = 360/alpha;
        if N_repeat_trial==floor(N_repeat_trial)
            N_repeat = N_repeat_trial;
        else
            N_repeat = round(N_repeat_trial);
            alpha = 360/N_repeat; %actual angle
            warning('Alpha angle was changed to %d',alpha);
            rho = R*(tand(alpha/2+45)^(1-NR));
        end
    end 
    MeshGeometry.Specimen.alpha = alpha; % store the relevant value for documentation.
%% Specimen geometry
        %% Create rectangular box with pointy edge
        sampleWidth = R;
        sampleThickness = 1;
        sampleHeight = H;
        numElementsWidth = NR;
        numElementsThickness = 1;
        if benchmark_flag
            numElementsThickness =N_repeat;
        end
        numElementsHeight = NH;
        %% Creating model geometry and mesh
        % A box is created with tri-linear hexahedral (hex8) elements using the
        % |hexMeshBox| function. The function offers the boundary faces with
        % seperate labels for the top, bottom, left, right, front, and back sides.
        % As such these can be used to define boundary conditions on the exterior. 

        % Create a box with hexahedral elements
        cubeDimensions=[sampleWidth sampleThickness sampleHeight]; %Dimensions
        cubeElementNumbers=[numElementsWidth numElementsThickness numElementsHeight]; %Number of elements
        outputStructType=2; %A structure compatible with mesh view
        [meshStruct]=hexMeshBox(cubeDimensions,cubeElementNumbers,outputStructType);

        %Access elements, nodes, and faces from the structure
        E=meshStruct.elements; %The elements 
        V=meshStruct.nodes; %The nodes (vertices)
        Fb=meshStruct.facesBoundary; %The boundary faces
        Cb=meshStruct.boundaryMarker; %The "colors" or labels for the boundary faces
        %% Algining specimen mesh (trnaslation, still cube)
        V(:,1) = V(:,1)+0.5*sampleWidth; 
        V(:,2) = V(:,2)+0.5*sampleThickness;
        V(:,3) = V(:,3)-0.5*sampleHeight;
        meshStruct.nodes = V;
        %% Applying mesh density bias
        x_layers = unique(V(:,1)); % rear layers by y val
        for i=1:length(x_layers)
            x_layer_nodelist{i} = find(V(:,1)==x_layers(i));
        end
        if ignore_formula
            V(:,1) = V(:,1).*(V(:,1)/R).^R_bias; % apply R bias
        else
            V(x_layer_nodelist{2},1) = rho;
            for i=3:numel(x_layer_nodelist)
                V(x_layer_nodelist{i},1) = V(x_layer_nodelist{i-1},1)*tand(alpha/2+45);
            end
        end
        V(:,3) = V(:,3).*abs(V(:,3)/H).^H_bias; % apply H bias
        %% Apply a rotation to parallel planes (layers in box) of nodes
        % identify the nodes in each layer parallel to x-z plane
        y_layers = unique(V(:,2)); % y values of each layer
        for i=1:length(y_layers)
            y_layer_nodelist{i} = find(V(:,2)==y_layers(i)); 
        end
        if (numel(y_layer_nodelist)-1)~=numElementsThickness
                    error('problem here');
        end
        rear_nodes = find(V(:,2)>0);
        % repeat specimen pattern for benchmark
        if benchmark_flag
            angle_increment = alpha;
        else 
            angle_increment = alpha/numElementsThickness;
        end
        for i=2:numel(y_layer_nodelist)
            % set x-y position of each layer to be a rotated copy (by alpha_increment degrees) of the previous layer
            V(y_layer_nodelist{i},1:2) = V(y_layer_nodelist{i-1},1:2)*[cosd(angle_increment), sind(angle_increment); -sind(angle_increment), cosd(angle_increment)]; 
        end
        %% convert hex8 to hex20 if needed, and partition mesh to hex and penta elements. 
        % E1 - hex elements
        % E2 - penta elements (created by merging nodes along the
        % degenerate facets of the hex, i.e. the centerline).
        if strcmp(elementType,'hex20')
            Fb_8 = Fb;
            Fb = [];
            [E,V,~,Fb_20]=hex8_hex20(E,V,{},Fb_8);
            Fb = Fb_20;
            meshStruct.elements=E;
            meshStruct.nodes=V;
            meshStruct.Fb=Fb;
            quadType='quad8';
            new_face_order = [1 3 5 7 2 4 6 8]; % compensate for GIBBON's quad8 nodal ordering scheme to allow work with FEBio
        else
            quadType='quad4';
            new_face_order = [1 2 3 4];
        end
        % merge centerline nodes
        [E,V,ind1,ind2]=mergeVertices(E,V);
        centerline_nodes = intersect(find(V(:,1)==0),(find(V(:,2)==0)));
        front_nodes = find(V(:,2)==0);
        % update Fb array with new indices.
        for i=1:size(Fb,1)
            for j=1:size(Fb,2)
%                 Cb(ind2(Fb(i,j))) = Cb(Fb(i,j));
                Fb(i,j) = ind2(Fb(i,j));
            end
        end
        % Partition elements E to hexaherdral (index '1') and pentaherdal (index '2')
        E1 = []; E1_ind = []; %hex8/hex20
        E2 = []; E2_ind = []; %penta6/penta15
        for i=1:size(E,1)
            if any(E(i,:)==centerline_nodes, 'all') % Element edge on centerline
                % this element should be converted to penta6
                E2(end+1,:)=E(i,:);
                E2_ind(end+1) = i;
            else % No edge centerline
                % this element should remain hex8/hex20
                E1(end+1,:)=E(i,:);
                E1_ind(end+1) = i;
            end
        end
        % Do the same for boundary faces and boundary faces colors
        Fb1_ind = []; Fb2_ind = []; Fb1 = [];
        Fb2 = []; Cb1 = []; Cb2 = [];
        for i=1:size(Fb,1)
            if any(Fb(i,:)==centerline_nodes, 'all')
                Fb2(end+1,:)=Fb(i,:);
                Fb2_ind(end+1) = i;
                Cb2(end+1,1) = Cb(Fb2_ind(end));
        %         Fb(Fb2_ind(end),4) = Fb(Fb2_ind(end),1);
        %         Fb(Fb2_ind(end),8) = Fb(Fb2_ind(end),5);
            else
                Fb1(end+1,:)=Fb(i,:);
                Fb1_ind(end+1) = i;
                Cb1(end+1,1) = Cb(Fb1_ind(end));
            end
        end
        % Dispose columns of nodes to convert hex to penta
        if strcmp(elementType,'hex20')
            % hex20 to penta15
            E2(:,[4,8,12,16,20]) = [];
        else
            % hex8 to penta6
            E2(:,[4,8]) = [];
        end
        % Move face 1 to interface between hexes and pentas
        Face1 = Fb(Cb==1,:);
        Fb2 = [Fb2;Face1+2]; 
        Cb2 = [Cb2;1*ones(size(Face1,1),1)];
        Fb1 = [Fb1;Face1+2]; 
        Cb1 = [Cb1;1*ones(size(Face1,1),1)];
        
        %% Update MeshGeometry structure
        MeshGeometry.Specimen.E = E; 
        MeshGeometry.Specimen.Fb = Fb;
        MeshGeometry.Specimen.Cb = Cb;
        MeshGeometry.Specimen.V = V;
        MeshGeometry.Specimen.nodes_id = 1:size(MeshGeometry.Specimen.V,1)';
        % Hexahedral
        MeshGeometry.Specimen.E_hex = E1; 
        MeshGeometry.Specimen.E_hex_ind = E1_ind; 
        MeshGeometry.Specimen.Fb_hex = Fb1;
        MeshGeometry.Specimen.Cb_hex = Cb1;
        % Pentahedral
        MeshGeometry.Specimen.E_penta = E2; 
        MeshGeometry.Specimen.E_penta_ind = E2_ind; 
        MeshGeometry.Specimen.Fb_penta = Fb2;
        MeshGeometry.Specimen.Cb_penta = Cb2;
        % 
        MeshGeometry.Specimen.numElementsThickness = numElementsThickness;
        % Specimen Facets lists
        MeshGeometry.Specimen.centerline_facets = MeshGeometry.Specimen.Fb(MeshGeometry.Specimen.Cb==1,new_face_order);
        MeshGeometry.Specimen.outer_facets = MeshGeometry.Specimen.Fb(MeshGeometry.Specimen.Cb==2,new_face_order);
        MeshGeometry.Specimen.front_facets = MeshGeometry.Specimen.Fb(MeshGeometry.Specimen.Cb==3,new_face_order);
        MeshGeometry.Specimen.back_facets = MeshGeometry.Specimen.Fb(MeshGeometry.Specimen.Cb==4,new_face_order);
        MeshGeometry.Specimen.bottom_facets = MeshGeometry.Specimen.Fb(MeshGeometry.Specimen.Cb==5,new_face_order);
        MeshGeometry.Specimen.top_facets = MeshGeometry.Specimen.Fb(MeshGeometry.Specimen.Cb==6,new_face_order);
        % Specimen nodes lists
        MeshGeometry.Specimen.top_surface_node_list = unique(MeshGeometry.Specimen.top_facets);        
        MeshGeometry.Specimen.centerline_nodes = centerline_nodes;
        MeshGeometry.Specimen.front_nodes = front_nodes;
        MeshGeometry.Specimen.bottom_nodes = unique(MeshGeometry.Specimen.bottom_facets);  
        MeshGeometry.Specimen.outer_nodes = unique(MeshGeometry.Specimen.outer_facets);
        MeshGeometry.Specimen.back_nodes = unique(MeshGeometry.Specimen.back_facets);
%% Indenter parameters
    contactInitialOffset =1e-3;
    numRefineStepsSphereIndenter = 5;
    IndenterHeight=1;  % height of the cylinder part
    pointSpacingIndenter=1;
    indenterRadius = Indenter.indenterRadius;
%     indentationMagnitude=20;
%     bcPrescribeMagnitudes=[0 0 -indentationMagnitude-contactInitialOffset]; %NB desired and effect of initial spacing
%     snr=30; % signal to noise ratio to add to the "experimental" data
    delta_alpha = 0; %0.2*alpha;
    alpha_lim = [0 alpha]+[-delta_alpha delta_alpha];
    phi_lim = [0 60];
    ns = 150;
    pointSpacing = x_layers(2)/5;
    resampleCurveOpt = 1;    
    interpMethod = 'linear';
%% Indenter geometry
        %%  Creating triangulated indenter surface model
        if benchmark_flag
            IndenterMeshInputStruct=struct;
            IndenterMeshInputStruct.sphereRadius=indenterRadius;
            IndenterMeshInputStruct.nRefine=numRefineStepsSphereIndenter;
            IndenterMeshInputStruct.cylinderHeight=IndenterHeight;
            IndenterMeshInputStruct.cylinderStepSize=pointSpacingIndenter;
            IndenterMeshInputStruct.patchType='tri';
            [E3,V3,~]=hemiSphereCylMesh(IndenterMeshInputStruct);
        else
            [E3,V3] = getIndenterMesh(alpha_lim, phi_lim, indenterRadius, ns, pointSpacing,resampleCurveOpt,interpMethod);
        end
        % offset indenter
        V3(:,3)=V3(:,3)+contactInitialOffset;
        center_of_mass=mean(V3,1);    
        %% combine specimen and indenter elements and nodes
        E3 = E3 + size(V,1); % increment indices
        E3_ind = (length(E2_ind)+length(E1_ind)+(1:1:size(E3,1))); 
        V = [V;V3]; % combined specimen and indenter nodal coordinates vector
        indenter_nodes_id = size(V,1)-size(V3,1)+1:size(V,1); % nodal IDs of indenter nodes -  V(indenter_nodes_id)==V3
        % Update MeshGeometry structure
        MeshGeometry.Indenter.V3 = V3;
        MeshGeometry.Indenter.E3 = E3;
        MeshGeometry.Indenter.E3_ind = E3_ind;
        MeshGeometry.Indenter.nodes_id = indenter_nodes_id;
        MeshGeometry.Indenter.center_of_mass = center_of_mass;
        MeshGeometry.Indenter.contactInitialOffset = contactInitialOffset;
        MeshGeometry.Indenter.radius = indenterRadius;
        MeshGeometry.Indenter.pointSpacingIndenter = pointSpacingIndenter;
        MeshGeometry.Indenter.IndenterHeight = IndenterHeight;
        MeshGeometry.Indenter.ns = ns;
        MeshGeometry.Indenter.pointSpacing = pointSpacing;
        MeshGeometry.Indenter.resampleCurveOpt = resampleCurveOpt;
        MeshGeometry.Indenter.interpMethod = interpMethod;
end
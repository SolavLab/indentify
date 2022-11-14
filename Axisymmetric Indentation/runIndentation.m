%% run_mesh_convergence_benchmark funciton
function [febio_spec,febioAnalysis, runFlag,savePath,MeshGeometry] = runIndentation(my_param,modelName,savePath)
    run_optimization = 0; % flag for carrying an identification study with synthetic test data (iFEA).
    FigStruct.Visible = 'off';
    output_all_nodes = 1;
    Specimen = my_param.Specimen;
    Indenter = my_param.Indenter;
    animationFlag=0;
    %% Specimen parameters
    MeshGeometry = axisymmetricGeometry(Specimen,Indenter);
    benchmark_flag = Specimen.benchmark_flag; % benchmark flag
    elementType = Specimen.elementType; % hex20 / hex8 
    pointSpacingGel = Specimen.R/Specimen.NR; % search radius factor for F2F and S2 contact formulations. currently not in use (only SE).
    %% Control parameters
    mkdir(savePath) % create 
    % Defining file names
    febioFebFileNamePart=sprintf('tempModel');
    febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
    febioLogFileName=[febioFebFileNamePart,'.txt']; %FEBio log file name
    febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
    febioLogFileName_stress=[febioFebFileNamePart,'_stress_out.txt']; %Log file name for exporting stress
    febioLogFileName_strain=[febioFebFileNamePart,'_strain_out.txt']; %Log file name for exporting strain
    febioLogFileName_pos=[febioFebFileNamePart,'_pos_out.txt']; %Log file name for exporting position of top surface nodes
    febioLogFileName_rigidBody=[febioFebFileNamePart,'_indenter_RB_out.txt'];
    %% Material parameter set
    % "true" Material parameter set
    p1_tru=my_param.p1; % 1st material parameter
    p2_tru=my_param.p2; % 2nd material parameter
    k_factor_tru=my_param.k_factor; %Bulk modulus factor
    %% FEA control settings
    nSteps=5; % number of indenations
    numTimeSteps_1=20; % time steps per indentation
    numTimeSteps=nSteps*numTimeSteps_1; %Number of time steps desired (3 indentations, 5 steps each)
    step_size=nSteps/numTimeSteps;
    max_refs=50; %Max reforms
    max_ups=10; %Set to zero to use full-Newton iterations
    opt_iter=10; %Optimum number of iterations
    max_retries=5; %Maximum number of retires
    dtmin=step_size/50; %Minimum time step size
    dtmax=step_size*2; %Maximum time step size
    %% Contact parametersP
    contactAlg=5;
    switch contactAlg
        case 1
            contactType='sticky';
        case 2
            contactType='facet-to-facet sliding';
        case 3
            contactType='sliding_with_gaps';
        case 4
            contactType='sliding2';
        case 5
            contactType ='sliding-elastic';
    end
    %% Plot settings
    fontSize1=15;
    fontSize2=25;
    faceAlpha1=0.8;
    faceAlpha2=0.3;
    markerSize1=20;
    markerSize2=6;
    lineWidth=3;
    %% Final positions of indenter
    random_depth = 0;
    nSteps=5; % excluding t=0;
    rotVec=zeros(nSteps,3);
    transVec=zeros(nSteps,3);
    %                 rotMat=zeros(3,3,nSteps);
            V3_steps=zeros([size(MeshGeometry.Indenter.V3) nSteps]);
    switch random_depth
        case 0 % predefined depth
            final_depth = -10;
            z_pos = linspace(0,final_depth,nSteps+1);
            transVec(:,3) = z_pos(2:end);
        case 2 % semi-random (same rng seed)
            rng(1,'twister'); %%%%<<<<<specify seed to generate *the same* random values for repeatability ~~!!!!!! 
            for ii=1:nSteps
                rotAxis=10*[rand-.5 rand-.5 rand+.5];
                rotAxis=rotAxis/norm(rotAxis);
                rotAng=0.1*rand;
                rotVec(ii,:)=rotAng*rotAxis;
%                     rotMat(:,:,ii)=rotationVectorToMatrix(rotVec(ii,:))';
                transVec(ii,:)=-[rand rand 2*ii+rand];
                V3_steps(:,:,ii)=MeshGeometry.Indenter.center_of_mass+(rotMat(:,:,ii)*(MeshGeometry.Indenter.V3-MeshGeometry.Indenter.center_of_mass)')'+transVec(ii,:);
            end
    end
    %% plot indentation positions
%     cFigure(FigStruct);
%     hold on;
%     gpatch(Fb,V,'y','y',faceAlpha2);
%     gpatch(E3,V3,'kw','none',faceAlpha2);
    colors=gjet(nSteps);
%     for ii=1:nSteps
%         gpatch(E3,V3_steps(:,:,ii),colors(ii,:),'none',faceAlpha2);
%     end
% 
%     axisGeom(gca,fontSize1);
%     camlight headlight
%     legendStrings={'cylinder','initial'};
%     for ii=1:nSteps
%         legendStrings=[legendStrings, {['indentation ', num2str(ii)]}];
%     end
%     legend(legendStrings);
        %% Plotting model boundary surfaces and a cut view
%     hFig=cFigure(FigStruct); hFig.Tag = 'temp';
% 
%     subplot(1,3,1); hold on; 
%     title('Model boundary surfaces and labels','FontSize',fontSize);
%     gpatch(Fb,V,Cb,'k',faceAlpha1); 
%     colormap(gjet(6)); icolorbar;
%     axisGeom(gca,fontSize);
% 
%     subplot(1,3,2); hold on; 
%     title('Model boundary surfaces and labels','FontSize',fontSize);
%     gpatch(Fb1,V,Cb1,'k',faceAlpha1); 
%     colormap(gjet(6)); icolorbar;
%     axisGeom(gca,fontSize);
% 
%     subplot(1,3,3); hold on; 
%     title('Model boundary surfaces and labels','FontSize',fontSize);
%     gpatch(Fb2,V,Cb2,'k',faceAlpha1); 
%     colormap(gjet(6)); icolorbar;
%     axisGeom(gca,fontSize);
%     hFig.Visible = 'off';
%     drawnow;
    %% Defining the boundary conditions
    % The visualization of the model boundary shows colors for each side of the
    % cube. These labels can be used to define boundary conditions. 

    bcSupportList= MeshGeometry.Specimen.bottom_nodes; 
    
    F_rear_wall=MeshGeometry.Specimen.back_facets; % Define rear wall constraint surfaces
    top_surface_facets = MeshGeometry.Specimen.top_facets;
    bcPrescribeList=MeshGeometry.Specimen.top_surface_node_list; %Node set part of selected face
    front_nodes = MeshGeometry.Specimen.front_nodes;
    centerline_nodes = MeshGeometry.Specimen.centerline_nodes;
    F_contact_master = MeshGeometry.Indenter.E3;
    indenter_nodes_id = MeshGeometry.Indenter.nodes_id;
        %% Visualizing boundary conditions. Markers plotted on the semi-transparent
    % model denote the nodes in the various boundary condition lists. 

%     hf=cFigure; hf.Tag = 'temp';
%     title('Boundary conditions','FontSize',fontSize);
%     xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
%     hold on;
% 
%     gpatch(Fb,V,'kw','k',0.5);
% 
%     hl(1)=plotV(V(bcSupportList,:),'k.','MarkerSize',markerSize);
%     hl(2)=plotV(V(bcPrescribeList,:),'r.','MarkerSize',markerSize);
% 
%     legend(hl,{'BC support','BC prescribe'});
% 
%     axisGeom(gca,fontSize);
%     camlight headlight; 
%     drawnow; 
    %% Defining the FEBio input structure
    % See also |febioStructTemplate| and |febioStruct2xml| and the FEBio user
    % manual.
        %% General
        %Get a template with default settings 
        [febio_spec]=febioStructTemplate;

        %febio_spec version 
        febio_spec.ATTR.version='3.0'; 

        %Module section
        febio_spec.Module.ATTR.type='solid'; 
        %% Control section
        febio_spec.Control = rmfield(febio_spec.Control, 'analysis');
        febio_spec.Control = rmfield(febio_spec.Control, 'solver');

        stepStruct.Control.time_steps=numTimeSteps;
        stepStruct.Control.step_size=step_size;
        stepStruct.Control.time_stepper.dtmin=dtmin;
        stepStruct.Control.time_stepper.max_retries=max_retries;
        stepStruct.Control.time_stepper.opt_iter=opt_iter;
        
        if isfield(stepStruct.Control.time_stepper, 'dtmax')
            stepStruct.Control.time_stepper = rmfield(stepStruct.Control.time_stepper, 'dtmax');
        end
        stepStruct.Control.plot_zero_state = 1;
        stepStruct.Control.output_level = 'OUTPUT_MUST_POINTS';
        stepStruct.Control.plot_level = 'PLOT_MUST_POINTS'; % PLOT_MUST_POINTS
        stepStruct.Control.time_stepper.dtmax.ATTR.lc=1;
        stepStruct.Control.time_stepper.dtmax.VAL =1;  
        stepStruct.Control.solver.symmetric_stiffness = 0;
        stepStruct.Control.solver.max_refs = max_refs;
        stepStruct.Control.solver.max_ups = max_ups;

        % febio_spec.Control.time_stepper.dtmax.ATTR.lc=1; % must point defined using loadcurve 1
        % febio_spec.Control.time_stepper.dtmax.VAL=dtmax; 

        %Add template based default settings to proposed control section
        [stepStruct.Control]=structComplete(stepStruct.Control,febio_spec.Control,1); %Complement provided with default if missing

        %Remove control field (part of template) since step specific control sections are used
        febio_spec=rmfield(febio_spec,'Control'); 

        %Step specific control section
        febio_spec.Step.step{1}.Control=stepStruct.Control;
        febio_spec.Step.step{1}.ATTR.id=1;
        %% Material section
        switch my_param.mat_type
            case 'Ogden_1st' % Ogden 1st order
                k_tru=0.5*p1_tru*k_factor_tru; %Bulk modulus = (initial shear modulus)x(k_factor)
                materialName1='Material1';
                febio_spec.Material.material{1}.ATTR.name=materialName1;
                febio_spec.Material.material{1}.ATTR.type='Ogden';
                febio_spec.Material.material{1}.ATTR.id=1;
                febio_spec.Material.material{1}.c1=p1_tru;
                febio_spec.Material.material{1}.m1=p2_tru;
                febio_spec.Material.material{1}.k=k_tru;
            case 'Ogden_symmetric' % Ogden 2nd order (symmetric)
                k_tru=p1_tru*k_factor_tru; %Bulk modulus = (initial shear modulus)x(k_factor)
                materialName1='Material1';
                febio_spec.Material.material{1}.ATTR.name=materialName1;
                febio_spec.Material.material{1}.ATTR.type='Ogden';
                febio_spec.Material.material{1}.ATTR.id=1;
                febio_spec.Material.material{1}.c1=p1_tru;
                febio_spec.Material.material{1}.m1=p2_tru;
                febio_spec.Material.material{1}.c2=p1_tru;
                febio_spec.Material.material{1}.m2=-p2_tru;
                febio_spec.Material.material{1}.k=k_tru;            
            case 'Mooney-Rivlin' % Mooney Rivlin         
                materialName1='Material1';
                k_tru=2*(p1_tru+p2_tru)*k_factor_tru; %Bulk modulus = (initial shear modulus)x(k_factor)
                febio_spec.Material.material{1}.ATTR.name=materialName1;
                febio_spec.Material.material{1}.ATTR.type='Mooney-Rivlin';
                febio_spec.Material.material{1}.ATTR.id=1;
                febio_spec.Material.material{1}.c1=p1_tru;
                febio_spec.Material.material{1}.c2=p2_tru;
%                 febio_spec.Material.material{1}.k=880*(febio_spec.Material.material{1}.c1+febio_spec.Material.material{1}.c2);
                febio_spec.Material.material{1}.k=k_tru;
            case 'Neo-Hookean(MR)' % Neo-Hookean (MR)
                materialName1='Material1';
                febio_spec.Material.material{1}.ATTR.name=materialName1;
                febio_spec.Material.material{1}.ATTR.type='Mooney-Rivlin';
                febio_spec.Material.material{1}.ATTR.id=1;
                febio_spec.Material.material{1}.c1=p1_tru;
                febio_spec.Material.material{1}.c2=0;
                febio_spec.Material.material{1}.k=(4/3)*p1_tru*(1-p2_tru)/(1-2*p2_tru);
            case 'Neo-Hookean' % Neo-Hookean 
                    materialName1='Material1';
                    febio_spec.Material.material{1}.ATTR.name=materialName1;
                    febio_spec.Material.material{1}.ATTR.type='neo-Hookean';
                    febio_spec.Material.material{1}.ATTR.id=1;
                    febio_spec.Material.material{1}.v=p2_tru;
                    febio_spec.Material.material{1}.E=2*p1_tru*(1+febio_spec.Material.material{1}.v);
            case 'Neo-Hookean(Young)' % Neo-Hookean (Young)
                    materialName1='Material1';
                    febio_spec.Material.material{1}.ATTR.name=materialName1;
                    febio_spec.Material.material{1}.ATTR.type='neo-Hookean';
                    febio_spec.Material.material{1}.ATTR.id=1;
                    febio_spec.Material.material{1}.v=p2_tru;
                    febio_spec.Material.material{1}.E=p1_tru;
        end
        % Indeter material (rigid)
        materialName2 = 'Material2';
        febio_spec.Material.material{2}.ATTR.name=materialName2;
        febio_spec.Material.material{2}.ATTR.type='rigid body';
        febio_spec.Material.material{2}.ATTR.id=2;
        febio_spec.Material.material{2}.density=1;
        febio_spec.Material.material{2}.center_of_mass=MeshGeometry.Indenter.center_of_mass;
        %% Mesh section
            % -> Nodes
                febio_spec.Mesh.Nodes{1}.ATTR.name='AllNodes'; %The node set name
                febio_spec.Mesh.Nodes{1}.node.ATTR.id=[MeshGeometry.Specimen.nodes_id,MeshGeometry.Indenter.nodes_id]'; %The node id's
                febio_spec.Mesh.Nodes{1}.node.VAL=[MeshGeometry.Specimen.V;MeshGeometry.Indenter.V3]; %The nodal coordinates
        if strcmp(elementType, 'hex20')
            % -> Elements
            partName1='Part1';
            febio_spec.Mesh.Elements{1}.ATTR.name=partName1; %Name of this part
            febio_spec.Mesh.Elements{1}.ATTR.type='hex20'; %Element type
            febio_spec.Mesh.Elements{1}.elem.ATTR.id=MeshGeometry.Specimen.E_hex_ind'; %Element id's
            febio_spec.Mesh.Elements{1}.elem.VAL=MeshGeometry.Specimen.E_hex; %The element matrix
            partName2='Part2';
            febio_spec.Mesh.Elements{2}.ATTR.name=partName2; %Name of this part
            febio_spec.Mesh.Elements{2}.ATTR.type='penta15'; %Element type
            febio_spec.Mesh.Elements{2}.elem.ATTR.id=MeshGeometry.Specimen.E_penta_ind'; %Element id's
            febio_spec.Mesh.Elements{2}.elem.VAL=MeshGeometry.Specimen.E_penta; %The element matrix

            partName3='Sphere';
            febio_spec.Mesh.Elements{3}.ATTR.name=partName3; %Name of this part
            febio_spec.Mesh.Elements{3}.ATTR.type='tri3'; %Element type of this set
            febio_spec.Mesh.Elements{3}.ATTR.mat=2; %material index for this set
            % febio_spec.Mesh.Elements{2}.ATTR.name='Sphere'; %Name of the element set
            febio_spec.Mesh.Elements{3}.elem.ATTR.id=MeshGeometry.Indenter.E3_ind'; %Element id's
            febio_spec.Mesh.Elements{3}.elem.VAL=MeshGeometry.Indenter.E3;

            % -> Surfaces
            surfaceName1='rear_wall_master';
            febio_spec.Mesh.Surface{1}.ATTR.name=surfaceName1;
            febio_spec.Mesh.Surface{1}.quad8.ATTR.id=(1:1:size(F_rear_wall,1))';
            febio_spec.Mesh.Surface{1}.quad8.VAL=F_rear_wall;

            surfaceName2='contact_slave_quad4';
            % top_surface_facets_tri3 = top_surface_facets(1,1:3);
            top_surface_facets_quad4 = top_surface_facets(1+MeshGeometry.Specimen.numElementsThickness:end,:); % actually is quad8
            included_indices = [];
            for i=1:size(top_surface_facets_quad4,1)
                if any(MeshGeometry.Specimen.V(top_surface_facets_quad4(i,:),1)<=MeshGeometry.Indenter.radius)
                    included_indices(end+1) = i;
                end
            end
            febio_spec.Mesh.Surface{2}.ATTR.name=surfaceName2;
            febio_spec.Mesh.Surface{2}.quad8.ATTR.id=(1:1:size(top_surface_facets_quad4(included_indices,:),1))'+febio_spec.Mesh.Surface{1}.quad8.ATTR.id(end);
            febio_spec.Mesh.Surface{2}.quad8.VAL=top_surface_facets_quad4(included_indices,:);
            % febio_spec.Mesh.Surface{2}.tri3.ATTR.id=(1:1:size(top_surface_facets_tri3,1))'+(size(F_rear_wall,1)+size(top_surface_facets_quad4,1));
            % febio_spec.Mesh.Surface{2}.tri3.VAL=top_surface_facets_tri3;


            surfaceName3='contact_master';
            febio_spec.Mesh.Surface{3}.ATTR.name=surfaceName3;
            febio_spec.Mesh.Surface{3}.tri3.ATTR.id=(1:1:size(F_contact_master,1))'+febio_spec.Mesh.Surface{2}.quad8.ATTR.id(end);
            febio_spec.Mesh.Surface{3}.tri3.VAL=F_contact_master;

            surfaceName4='contact_slave_tri3';
            % top_surface_facets_tri3 = top_surface_facets(1,1:3);
            top_surface_facets_tri3 = top_surface_facets(1:MeshGeometry.Specimen.numElementsThickness,[1:3,5:7]);
            febio_spec.Mesh.Surface{4}.ATTR.name=surfaceName4;
            febio_spec.Mesh.Surface{4}.tri6.ATTR.id=(1:1:size(top_surface_facets_tri3,1))'+febio_spec.Mesh.Surface{3}.tri3.ATTR.id(end);
            febio_spec.Mesh.Surface{4}.tri6.VAL=top_surface_facets_tri3;
        else 
             % -> Elements
            partName1='Part1';
            febio_spec.Mesh.Elements{1}.ATTR.name=partName1; %Name of this part
            febio_spec.Mesh.Elements{1}.ATTR.type='hex8'; %Element type
            febio_spec.Mesh.Elements{1}.elem.ATTR.id=MeshGeometry.Specimen.E_hex_ind'; %Element id's
            febio_spec.Mesh.Elements{1}.elem.VAL=MeshGeometry.Specimen.E_hex; %The element matrix
            partName2='Part2';
            febio_spec.Mesh.Elements{2}.ATTR.name=partName2; %Name of this part
            febio_spec.Mesh.Elements{2}.ATTR.type='penta6'; %Element type
            febio_spec.Mesh.Elements{2}.elem.ATTR.id=MeshGeometry.Specimen.E_penta_ind'; %Element id's
            febio_spec.Mesh.Elements{2}.elem.VAL=MeshGeometry.Specimen.E_penta; %The element matrix

            partName3='Sphere';
            febio_spec.Mesh.Elements{3}.ATTR.name=partName3; %Name of this part
            febio_spec.Mesh.Elements{3}.ATTR.type='tri3'; %Element type of this set
            febio_spec.Mesh.Elements{3}.ATTR.mat=2; %material index for this set
            % febio_spec.Mesh.Elements{2}.ATTR.name='Sphere'; %Name of the element set
            febio_spec.Mesh.Elements{3}.elem.ATTR.id=MeshGeometry.Indenter.E3_ind'; %Element id's
            febio_spec.Mesh.Elements{3}.elem.VAL=MeshGeometry.Indenter.E3;
            

            % -> Surfaces
            surfaceName1='rear_wall_master';
            febio_spec.Mesh.Surface{1}.ATTR.name=surfaceName1;
            febio_spec.Mesh.Surface{1}.quad4.ATTR.id=(1:1:size(F_rear_wall,1))';
            febio_spec.Mesh.Surface{1}.quad4.VAL=F_rear_wall;

            surfaceName2='contact_slave_quad4';
            % top_surface_facets_tri3 = top_surface_facets(1,1:3);
            top_surface_facets_quad4 = top_surface_facets(1+MeshGeometry.Specimen.numElementsThickness:end,:);
            included_indices = [];
            for i=1:size(top_surface_facets_quad4,1)
                if any(MeshGeometry.Specimen.V(top_surface_facets_quad4(i,:),1)<=MeshGeometry.Indenter.radius)
                    included_indices(end+1) = i;
                end
            end
            febio_spec.Mesh.Surface{2}.ATTR.name=surfaceName2;
            febio_spec.Mesh.Surface{2}.quad4.ATTR.id=(1:1:size(top_surface_facets_quad4(included_indices,:),1))'+febio_spec.Mesh.Surface{1}.quad4.ATTR.id(end);
            febio_spec.Mesh.Surface{2}.quad4.VAL=top_surface_facets_quad4(included_indices,:);
            % febio_spec.Mesh.Surface{2}.tri3.ATTR.id=(1:1:size(top_surface_facets_tri3,1))'+(size(F_rear_wall,1)+size(top_surface_facets_quad4,1));
            % febio_spec.Mesh.Surface{2}.tri3.VAL=top_surface_facets_tri3;


            surfaceName3='contact_master';
            febio_spec.Mesh.Surface{3}.ATTR.name=surfaceName3;
            febio_spec.Mesh.Surface{3}.tri3.ATTR.id=(1:1:size(F_contact_master,1))'+febio_spec.Mesh.Surface{2}.quad4.ATTR.id(end);
            febio_spec.Mesh.Surface{3}.tri3.VAL=F_contact_master;

            surfaceName4='contact_slave_tri3';
            % top_surface_facets_tri3 = top_surface_facets(1,1:3);
            top_surface_facets_tri3 = top_surface_facets(1:MeshGeometry.Specimen.numElementsThickness,1:3);
            febio_spec.Mesh.Surface{4}.ATTR.name=surfaceName4;
            febio_spec.Mesh.Surface{4}.tri3.ATTR.id=(1:1:size(top_surface_facets_tri3,1))'+febio_spec.Mesh.Surface{3}.tri3.ATTR.id(end);
            febio_spec.Mesh.Surface{4}.tri3.VAL=top_surface_facets_tri3;
        end
        % -> Surface pairs
        contactPairName{1}='Contact1';
        febio_spec.Mesh.SurfacePair{1}.ATTR.name=contactPairName{1};
        febio_spec.Mesh.SurfacePair{1}.primary=surfaceName2; 
        febio_spec.Mesh.SurfacePair{1}.secondary=surfaceName3;
        
        contactPairName{2}='Contact2';
        febio_spec.Mesh.SurfacePair{2}.ATTR.name=contactPairName{2};
        febio_spec.Mesh.SurfacePair{2}.primary=surfaceName4; 
        febio_spec.Mesh.SurfacePair{2}.secondary=surfaceName3;
        % -> NodeSets
        nodeSetName1='bcSupportList';
        nodeSetName2='bcPrescribeList';
        nodeSetName3='bcFrontSymmetryList';
        nodeSetName4 = 'BcCenterline';
        nodeSetName5 = 'indenterRigid';

        febio_spec.Mesh.NodeSet{1}.ATTR.name=nodeSetName1;
        febio_spec.Mesh.NodeSet{1}.node.ATTR.id=bcSupportList(:);

        febio_spec.Mesh.NodeSet{2}.ATTR.name=nodeSetName2;
        febio_spec.Mesh.NodeSet{2}.node.ATTR.id=bcPrescribeList(:);

        febio_spec.Mesh.NodeSet{3}.ATTR.name=nodeSetName3;
        febio_spec.Mesh.NodeSet{3}.node.ATTR.id=front_nodes(:);

        febio_spec.Mesh.NodeSet{4}.ATTR.name=nodeSetName4;
        febio_spec.Mesh.NodeSet{4}.node.ATTR.id=centerline_nodes(:);

        febio_spec.Mesh.NodeSet{5}.ATTR.name=nodeSetName5;
        febio_spec.Mesh.NodeSet{5}.node.ATTR.id=indenter_nodes_id(:);
        %% MeshDomains section
        febio_spec.MeshDomains.SolidDomain{1}.ATTR.name=partName1;
        febio_spec.MeshDomains.SolidDomain{1}.ATTR.mat=materialName1;

        febio_spec.MeshDomains.SolidDomain{2}.ATTR.name=partName2;
        febio_spec.MeshDomains.SolidDomain{2}.ATTR.mat=materialName1;

        febio_spec.MeshDomains.ShellDomain.ATTR.name=partName3;
        febio_spec.MeshDomains.ShellDomain.ATTR.mat=materialName2;
        %% Boundary condition section 
        %-> Fix boundary conditions
        febio_spec.Boundary.bc{1}.ATTR.type='fix';
        febio_spec.Boundary.bc{1}.ATTR.node_set=nodeSetName1;
        febio_spec.Boundary.bc{1}.dofs='x,y,z';
        if ~benchmark_flag          
            febio_spec.Boundary.bc{2}.ATTR.type='fix';
            febio_spec.Boundary.bc{2}.ATTR.node_set=nodeSetName4;
            febio_spec.Boundary.bc{2}.dofs='x,y';
            febio_spec.Boundary.bc{3}.ATTR.type='fix';
            febio_spec.Boundary.bc{3}.ATTR.node_set=nodeSetName3;
            febio_spec.Boundary.bc{3}.dofs='y';
        end
        %% Output section 
        % -> log file
        % nodal displacements
        febio_spec.Output.logfile.ATTR.file=febioLogFileName;
        febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
        febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
        febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';
            %  febio_spec.Output.logfile.node_data{1}.VAL=1:size(V,1); % all nodes
            %  febio_spec.Output.logfile.node_data{1}.VAL=top_surface_node_list'; % Only top nodes
            %  febio_spec.Output.logfile.node_data{1}.VAL =1:size(MeshGeometry.Specimen.V,1); % only specimen nodes
        %  nodal coordinates
        febio_spec.Output.logfile.node_data{2}.ATTR.file=febioLogFileName_pos;
        febio_spec.Output.logfile.node_data{2}.ATTR.data='x;y;z';
        febio_spec.Output.logfile.node_data{2}.ATTR.delim=',';
            %  febio_spec.Output.logfile.node_data{2}.VAL=top_surface_node_list';
            %  febio_spec.Output.logfile.node_data{2}.VAL=1:size(MeshGeometry.Specimen.V,1); % only specimen nodes
        %  rigid body data (foreces, moments and force)
        febio_spec.Output.logfile.rigid_body_data{1}.ATTR.file=febioLogFileName_rigidBody;
        febio_spec.Output.logfile.rigid_body_data{1}.ATTR.data='Fx;Fy;Fz;My;Mz;Mz;z';
        febio_spec.Output.logfile.rigid_body_data{1}.ATTR.delim=',';
        % Elemental stresses
        febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_stress;
        febio_spec.Output.logfile.element_data{1}.ATTR.data='s1;s2;s3';
        febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';
        febio_spec.Output.logfile.element_data{1}.VAL=1:size(MeshGeometry.Specimen.E,1);
        % Elemental strains
        febio_spec.Output.logfile.element_data{2}.ATTR.file=febioLogFileName_strain;
        febio_spec.Output.logfile.element_data{2}.ATTR.data='E1;E2;E3;Exy;J';
        febio_spec.Output.logfile.element_data{2}.ATTR.delim=',';
        febio_spec.Output.logfile.element_data{2}.VAL=1:size(MeshGeometry.Specimen.E,1);
        %% Constraints
        if ~benchmark_flag
            febio_spec.Constraints.constraint{1}.ATTR.type = 'symmetry plane';
            febio_spec.Constraints.constraint{1}.ATTR.surface = surfaceName1;
            febio_spec.Constraints.constraint{1}.laugon = 0; %%%%%%%% Setting this to 1 increased runtime from 10.5 to 23 seconds
            febio_spec.Constraints.constraint{1}.penalty = 1e6;
            febio_spec.Constraints.constraint{1}.tol = 1e-8;
            febio_spec.Constraints.constraint{1}.maxaug = 50;
        end
        %% Contact section 
        % Contact between indentor and specimen
        % Currently involves 2 identical contact cards (1 extra exclusivly for the
        % tri3 surface for the tri3/tri6 at r=0)
        contact_penalty = 5;
        auto_penalty = 1;
        two_pass_flag = 0;
        laugon = 1;
        for i=1:numel(febio_spec.Mesh.SurfacePair)
            switch contactType
                case 'sticky'
                    febio_spec.Contact.contact{i}.ATTR.surface_pair=contactPairName{i};
                    febio_spec.Contact.contact{i}.ATTR.type='sticky';
                    febio_spec.Contact.contact{i}.penalty=10;
                    febio_spec.Contact.contact{i}.laugon=1;
                    febio_spec.Contact.contact{i}.tolerance=0.1;
                    febio_spec.Contact.contact{i}.minaug=0;
                    febio_spec.Contact.contact{i}.maxaug=10;
                    febio_spec.Contact.contact{i}.snap_tol=0;
                    febio_spec.Contact.contact{i}.max_traction=0;
                    febio_spec.Contact.contact{i}.search_tolerance=0.1;
                case 'facet-to-facet sliding'
                    febio_spec.Contact.contact{i}.ATTR.surface_pair=contactPairName{i};
                    febio_spec.Contact.contact{i}.ATTR.type='facet-to-facet sliding';
                    febio_spec.Contact.contact{i}.penalty=200;
                    febio_spec.Contact.contact{i}.auto_penalty=1;
                    febio_spec.Contact.contact{i}.two_pass=0;
                    febio_spec.Contact.contact{i}.laugon=0;
                    febio_spec.Contact.contact{i}.tolerance=0.1;
                    febio_spec.Contact.contact{i}.gaptol=0;
                    febio_spec.Contact.contact{i}.minaug=0;
                    febio_spec.Contact.contact{i}.maxaug=10;
                    febio_spec.Contact.contact{i}.search_tol=0.01;
                    febio_spec.Contact.contact{i}.search_radius=mean(pointSpacingGel)/2;
                case 'sliding_with_gaps'
                    febio_spec.Contact.contact{i}.ATTR.surface_pair=contactPairName{i};
                    febio_spec.Contact.contact{i}.ATTR.type='sliding_with_gaps';
                    febio_spec.Contact.contact{i}.penalty=100;
                    febio_spec.Contact.contact{i}.auto_penalty=1;
                    febio_spec.Contact.contact{i}.two_pass=0;
                    febio_spec.Contact.contact{i}.laugon=0;
                    febio_spec.Contact.contact{i}.tolerance=0.1;
                    febio_spec.Contact.contact{i}.gaptol=0;
                    febio_spec.Contact.contact{i}.minaug=0;
                    febio_spec.Contact.contact{i}.maxaug=10;
                    febio_spec.Contact.contact{i}.fric_coeff=0;
                    febio_spec.Contact.contact{i}.fric_penalty=0;
                    febio_spec.Contact.contact{i}.ktmult=1;
                    febio_spec.Contact.contact{i}.seg_up=0;
                    febio_spec.Contact.contact{i}.search_tol=0.01;
                case 'sliding2'
                    febio_spec.Contact.contact{i}.ATTR.surface_pair=contactPairName{i};
                    febio_spec.Contact.contact{i}.ATTR.type='sliding2';
                    febio_spec.Contact.contact{i}.penalty=30;
                    febio_spec.Contact.contact{i}.auto_penalty=1;
                    febio_spec.Contact.contact{i}.two_pass=0;
                    febio_spec.Contact.contact{i}.laugon=0;
                    febio_spec.Contact.contact{i}.tolerance=0.1;
                    febio_spec.Contact.contact{i}.gaptol=0;
                    febio_spec.Contact.contact{i}.symmetric_stiffness=0;
                    febio_spec.Contact.contact{i}.search_tol=0.01;
                    febio_spec.Contact.contact{i}.search_radius=mean(pointSpacingGel)/2;
                case 'sliding-elastic'
                    febio_spec.Contact.contact{i}.ATTR.surface_pair=contactPairName{i};
                    febio_spec.Contact.contact{i}.ATTR.type='sliding-elastic';
                    febio_spec.Contact.contact{i}.penalty=contact_penalty;
                    febio_spec.Contact.contact{i}.auto_penalty=auto_penalty;
                    febio_spec.Contact.contact{i}.two_pass=two_pass_flag;
                    febio_spec.Contact.contact{i}.laugon=laugon;
                    febio_spec.Contact.contact{i}.tolerance=0.01;
                    febio_spec.Contact.contact{i}.minaug=0;
                    febio_spec.Contact.contact{i}.maxaug=10;
                    febio_spec.Contact.contact{i}.search_tol=0.1;
                    febio_spec.Contact.contact{i}.search_radius=1;
                    febio_spec.Contact.contact{i}.fric_coeff =1E8;
                    febio_spec.Contact.contact{i}.symmetric_stiffness=0;
            end
        end
        %% Step
        ii=1;
        febio_spec.Step.step{ii}.Rigid.rigid_constraint{1}.ATTR.name='RigidPrescribeX';
        febio_spec.Step.step{ii}.Rigid.rigid_constraint{1}.ATTR.type='prescribe';
        febio_spec.Step.step{ii}.Rigid.rigid_constraint{1}.rb=2;
        febio_spec.Step.step{ii}.Rigid.rigid_constraint{1}.dof='Rx';
        febio_spec.Step.step{ii}.Rigid.rigid_constraint{1}.value.ATTR.lc=2;
        febio_spec.Step.step{ii}.Rigid.rigid_constraint{1}.value.VAL=0; %1
        febio_spec.Step.step{ii}.Rigid.rigid_constraint{1}.relative=0;


        febio_spec.Step.step{ii}.Rigid.rigid_constraint{2}.ATTR.name='RigidPrescribeY';
        febio_spec.Step.step{ii}.Rigid.rigid_constraint{2}.ATTR.type='prescribe';
        febio_spec.Step.step{ii}.Rigid.rigid_constraint{2}.rb=2;
        febio_spec.Step.step{ii}.Rigid.rigid_constraint{2}.dof='Ry';
        febio_spec.Step.step{ii}.Rigid.rigid_constraint{2}.value.ATTR.lc=3;
        febio_spec.Step.step{ii}.Rigid.rigid_constraint{2}.value.VAL=0; %1
        febio_spec.Step.step{ii}.Rigid.rigid_constraint{2}.relative=0;

        febio_spec.Step.step{ii}.Rigid.rigid_constraint{3}.ATTR.name='RigidPrescribeZ';
        febio_spec.Step.step{ii}.Rigid.rigid_constraint{3}.ATTR.type='prescribe';
        febio_spec.Step.step{ii}.Rigid.rigid_constraint{3}.rb=2;
        febio_spec.Step.step{ii}.Rigid.rigid_constraint{3}.dof='Rz';
        febio_spec.Step.step{ii}.Rigid.rigid_constraint{3}.value.ATTR.lc=4;
        febio_spec.Step.step{ii}.Rigid.rigid_constraint{3}.value.VAL=1;
        febio_spec.Step.step{ii}.Rigid.rigid_constraint{3}.relative=0;

        febio_spec.Step.step{ii}.Rigid.rigid_constraint{4}.ATTR.name='RigidPrescribeRx';
        febio_spec.Step.step{ii}.Rigid.rigid_constraint{4}.ATTR.type='prescribe';
        febio_spec.Step.step{ii}.Rigid.rigid_constraint{4}.rb=2;
        febio_spec.Step.step{ii}.Rigid.rigid_constraint{4}.dof='Ru';
        febio_spec.Step.step{ii}.Rigid.rigid_constraint{4}.value.ATTR.lc=5;
        febio_spec.Step.step{ii}.Rigid.rigid_constraint{4}.value.VAL=0; %1
        febio_spec.Step.step{ii}.Rigid.rigid_constraint{4}.relative=0;

        febio_spec.Step.step{ii}.Rigid.rigid_constraint{5}.ATTR.name='RigidPrescribeRy';
        febio_spec.Step.step{ii}.Rigid.rigid_constraint{5}.ATTR.type='prescribe';
        febio_spec.Step.step{ii}.Rigid.rigid_constraint{5}.rb=2;
        febio_spec.Step.step{ii}.Rigid.rigid_constraint{5}.dof='Rv';
        febio_spec.Step.step{ii}.Rigid.rigid_constraint{5}.value.ATTR.lc=6;
        febio_spec.Step.step{ii}.Rigid.rigid_constraint{5}.value.VAL=0; %1
        febio_spec.Step.step{ii}.Rigid.rigid_constraint{5}.relative=0;

        febio_spec.Step.step{ii}.Rigid.rigid_constraint{6}.ATTR.name='RigidPrescribeRz';
        febio_spec.Step.step{ii}.Rigid.rigid_constraint{6}.ATTR.type='prescribe';
        febio_spec.Step.step{ii}.Rigid.rigid_constraint{6}.rb=2;
        febio_spec.Step.step{ii}.Rigid.rigid_constraint{6}.dof='Rw';
        febio_spec.Step.step{ii}.Rigid.rigid_constraint{6}.value.ATTR.lc=7;
        febio_spec.Step.step{ii}.Rigid.rigid_constraint{6}.value.VAL=0; %1
        febio_spec.Step.step{ii}.Rigid.rigid_constraint{6}.relative=0;
        %% load controllers
        febio_spec.LoadData.load_controller{1}.ATTR.id=1;
        febio_spec.LoadData.load_controller{1}.ATTR.type='loadcurve';
        febio_spec.LoadData.load_controller{1}.interpolate='STEP';
        ordinate_data = dtmax*ones(nSteps+1,1); ordinate_data(1) = 0;
        febio_spec.LoadData.load_controller{1}.points.point.VAL=[(0:nSteps)' ordinate_data];

        % load translations and rotations 
        febio_spec.LoadData.load_controller{2}.ATTR.id=2;
        febio_spec.LoadData.load_controller{2}.ATTR.type='loadcurve';
        febio_spec.LoadData.load_controller{2}.interpolate='LINEAR';
        febio_spec.LoadData.load_controller{2}.points.point.VAL=[(0:nSteps)' [0; transVec(:,1)]];

        febio_spec.LoadData.load_controller{3}.ATTR.id=3;
        febio_spec.LoadData.load_controller{3}.ATTR.type='loadcurve';
        febio_spec.LoadData.load_controller{3}.interpolate='LINEAR';
        febio_spec.LoadData.load_controller{3}.points.point.VAL=[(0:nSteps)' [0; transVec(:,2)]];

        febio_spec.LoadData.load_controller{4}.ATTR.id=4;
        febio_spec.LoadData.load_controller{4}.ATTR.type='loadcurve';
        febio_spec.LoadData.load_controller{4}.interpolate='LINEAR';
        febio_spec.LoadData.load_controller{4}.points.point.VAL=[(0:nSteps)' [0; transVec(:,3)]];

        febio_spec.LoadData.load_controller{5}.ATTR.id=5;
        febio_spec.LoadData.load_controller{5}.ATTR.type='loadcurve';
        febio_spec.LoadData.load_controller{5}.interpolate='LINEAR';
        febio_spec.LoadData.load_controller{5}.points.point.VAL=[(0:nSteps)' [0; rotVec(:,1)]];

        febio_spec.LoadData.load_controller{6}.ATTR.id=6;
        febio_spec.LoadData.load_controller{6}.ATTR.type='loadcurve';
        febio_spec.LoadData.load_controller{6}.interpolate='LINEAR';
        febio_spec.LoadData.load_controller{6}.points.point.VAL=[(0:nSteps)' [0; rotVec(:,2)]];

        febio_spec.LoadData.load_controller{7}.ATTR.id=7;
        febio_spec.LoadData.load_controller{7}.ATTR.type='loadcurve';
        febio_spec.LoadData.load_controller{7}.interpolate='LINEAR';
        febio_spec.LoadData.load_controller{7}.points.point.VAL=[(0:nSteps)' [0; rotVec(:,3)]];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Exporting the FEBio input file
    % Exporting the febio_spec structure to an FEBio input file is done using
    % the |febioStruct2xml| function. 
    run_attempts = 5;
    runFlag = 0;
    for ii=1:run_attempts
        febioStruct2xml(febio_spec,febioFebFileName); %Exporting to file and domNode
        %% Running the FEBio analysis
        % To run the analysis defined by the created FEBio input file the
        % |runMonitorFEBio| function is used. The input for this function is a
        % structure defining job settings e.g. the FEBio input file name. The
        % optional output runFlag informs the user if the analysis was run
        % succesfully. 

        febioAnalysis.run_filename=febioFebFileName; %The input file name
        febioAnalysis.run_logname=febioLogFileName; %The name for the log file
        febioAnalysis.disp_on=1; %Display information on the command window
        % febioAnalysis.runMode='internal';
        febioAnalysis.runMode=my_param.runMode;
        febioAnalysis.maxLogCheckTime=300; %Max log file checking time
        sprintf('Executing attempt number %d/%d', ii,run_attempts)
        [runFlag]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!    
        if runFlag % Successful job 
            sprintf('Attempt number %d successful', ii);
            if animationFlag
                %% Import FEBio results 
                % Importing nodal displacements from a log file
                [time_mat, N_disp_mat,~]=importFEBio_logfile(fullfile(savePath,febioLogFileName_disp)); %Nodal displacements
                time_mat=[0; time_mat(:)]; %Time

                N_disp_mat=N_disp_mat(:,2:end,:);
                sizImport=size(N_disp_mat);
                sizImport(3)=sizImport(3)+1;
                N_disp_mat_n=zeros(sizImport);
                N_disp_mat_n(:,:,2:end)=N_disp_mat;
                N_disp_mat=N_disp_mat_n;
                N_disp_magnitude_mat=squeeze(sqrt(sum(N_disp_mat.^2,2)));
                DN=N_disp_mat(:,:,end);
                DN_magnitude=N_disp_magnitude_mat(:,end);
                V = [MeshGeometry.Specimen.V;MeshGeometry.Indenter.V3];
                V_def=V+DN;
                V_DEF=N_disp_mat+repmat(V,[1 1 size(N_disp_mat,3)]);
                X_DEF=V_DEF(:,1,:);
                Y_DEF=V_DEF(:,2,:);
                Z_DEF=V_DEF(:,3,:);
                Fb1 = MeshGeometry.Specimen.Fb;
                center_of_mass = MeshGeometry.Indenter.center_of_mass;
%                 E2_joined = E3;
                [CF]=vertexToFaceMeasure(Fb1,DN_magnitude);
                % Importing rigid body data from a log file
                [~, RigidBody_mat,~]=importFEBio_logfile(fullfile(savePath,febioLogFileName_rigidBody)); %Nodal displacements
                RigidBody_mat = squeeze(RigidBody_mat)';
                RigidBody_mat(:,1) = []; %remove Indices
                RigidBody_mat = [zeros(1,size(RigidBody_mat,2)); RigidBody_mat]; % Add zero state row
                RigidBody_displacement=abs(RigidBody_mat(:,3)-center_of_mass(3));
                RigidBody_force_magnitude=sqrt(sum(RigidBody_mat(:,1:3).^2,2));
                RigidBody_Z_force_magnitude=abs(RigidBody_mat(:,3));
                RigidBody_torque_magnitude=sqrt(sum(RigidBody_mat(:,4:6).^2,2));
                 % Importing surface nodes nodal position from a log file
                [time_mat, N_pos_mat,~]=importFEBio_logfile(fullfile(savePath,febioLogFileName_pos),1,1); %Nodal position
                N_pos_mat(:,:,1) = V; % Get zero state from febio_spec
                %% Plotting the simulated results using |anim8| to visualize and animate
                % deformations
                colors=gjet(nSteps);
                % Create basic view and store graphics handle to initiate animation
                hf=cFigure; %Open figure
                hp1=gpatch(Fb1,V_def,CF,'none',1); %Add graphics object to animate
                % fake legend
                for ii=1:nSteps
                        gpatch([1 2 3],zeros(3),colors(ii,:),'none',faceAlpha2); %Add graphics object to animate
                end
                % title with force result
                ht1=gtitle([num2str(RigidBody_Z_force_magnitude(end),'%.2f') ' N/',num2str(MeshGeometry.Specimen.alpha),char(hex2dec('2da'))],fontSize2);
                % plot indenter
                E3 = MeshGeometry.Indenter.E3;
                hp2=gpatch(E3,V_def,colors(ii,:),'none',faceAlpha2); %Add graphics object to animate

                axisGeom(gca,fontSize1);
                colormap cool;
                hc=colorbar;
                title(hc,'dispMgn [mm]');
                caxis([0 max(N_disp_magnitude_mat(:))]);
                axis([min(X_DEF(:)) max(X_DEF(:)) min(Y_DEF(:)) max(Y_DEF(:)) min(Z_DEF(:)) max(Z_DEF(:))]);
                camlight headlight;

                % Set up animation features
                animStruct.Time=time_mat; %The time vector
                for qt=1:1:size(N_disp_mat,3) %Loop over time increments
                    DN=N_disp_mat(:,:,qt); %Current displacement
                    DN_magnitude=N_disp_magnitude_mat(:,qt); %Current displacement magnitude
                    V_def=V+DN; %Current nodal coordinates
                    [CF]=vertexToFaceMeasure(Fb1,DN_magnitude); %Current color data to use
                    % indenter color according to time step
                    colorInd=ceil(time_mat(qt));
                    colorInd(colorInd==0)=1;

                    %Set entries in animation structure
                    animStruct.Handles{qt}=[hp1 hp1 hp2 hp2 ht1]; %Handles of objects to animate
                    animStruct.Props{qt}={'Vertices','CData','Vertices','FaceColor','String'}; %Properties of objects to animate
                    animStruct.Set{qt}={V_def,CF,V_def,colors(colorInd,:),[num2str(RigidBody_Z_force_magnitude(qt),'%.2f') ' N/',num2str(MeshGeometry.Specimen.alpha),char(hex2dec('2da'))]}; %Property values for to set in order to animate
                end
                anim8(hf,animStruct); %Initiate animation feature
                set(gca,'FontSize',fontSize2);
                legendStrings={'cylinder'};
                for ii=1:nSteps
                    legendStrings=[legendStrings, {['indentation ', num2str(ii)]}];
                end
                legend(legendStrings);
            %     addFigureButtons
                drawnow;
            end
            break;
        else % Failed job - reduce contact penalty by 20% and retry
            sprintf('Attempt number %d failed', ii);
            for kk=1:numel(febio_spec.Contact.contact{i})
                febio_spec.Contact.contact{kk}.penalty = febio_spec.Contact.contact{kk}.penalty*0.8;
            end
        end
    end
    %% Close temporary figures
    close(findobj('Type', 'figure','Tag','temp'));
    %% run parameter identification problem <<< Moved to seperate script
end
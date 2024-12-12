%% NOTICE: the following code is an adaptation of DEMO_febio_0006_sphere_indentation
%Copyright (C) 2006-2023 Kevin Mattheus Moerman and the GIBBON contributors,
%taken from the GIBBON Toolbox (www.gibboncode.org) under the license
%provided therein:
%(https://github.com/gibbonCode/GIBBON/blob/master/LICENSE).

%This function runs an indentation simulation and saves the following data at set time steps alone:
% Indentor force
% Displacement and position of all nodes at the top surface


%% Setup axyisymmetric indentation FE model and run with FEBio.
function [febio_spec,febioAnalysis,runFlag] = runAnisotropicIndentation(my_param,disp_on)

%% Control parameters
savePath = my_param.savePath;
if ~exist(savePath,"dir")
    mkdir(savePath) % create
end
% Defining file names
febioFebFileNamePart='tempModel';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=[febioFebFileNamePart,'.txt']; %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_pos=[febioFebFileNamePart,'_pos_out.txt']; %Log file name for exporting position of top surface nodes
febioLogFileName_rigidBody=[febioFebFileNamePart,'_indenter_RB_out.txt'];

%% Specimen parameters
%Access elements, nodes, and faces from the structure
MeshGeometry=my_param.MeshGeometry;
E1=MeshGeometry.Specimen.elements; %The elements
V1=MeshGeometry.Specimen.nodes; %The nodes (vertices)
Fb1=MeshGeometry.Specimen.facesBoundary; %The boundary faces
Cb1=MeshGeometry.Specimen.boundaryMarker; %The "colors" or labels for the boundary faces
E2=MeshGeometry.Indenter.elements; %The elements
V2=MeshGeometry.Indenter.nodes; %The nodes (vertices)
elementType = MeshGeometry.Specimen.elementType; % hex20 / hex8

%Define applied displacement
sphereDisplacement=my_param.sphereDisplacement;

%% Material parameter set
%Each material type will have a different amount of parameters to consider.
%It is assumed the anisotropy is always in the same direction [0,1,0].
switch my_param.mat_type
    case 'trans iso Mooney-Rivlin'
        c1 = my_param.matParameters(1);
        c2 = my_param.matParameters(2);
        c3 = my_param.matParameters(3);
        c4 = my_param.matParameters(4);
        c5 = my_param.matParameters(5);
        lam_max = my_param.matParameters(6);
        k = my_param.matParameters(7);
    case 'trans iso Veronda-Westmann'
        c1 = my_param.matParameters(1);
        c2 = my_param.matParameters(2);
        c3 = my_param.matParameters(3);
        c4 = my_param.matParameters(4);
        c5 = my_param.matParameters(5);
        lam_max = my_param.matParameters(6);
        k = my_param.matParameters(7);
    case 'muscle material'
        g1 = my_param.matParameters(1);
        g2 = my_param.matParameters(2);
        p1 = my_param.matParameters(3);
        p2 = my_param.matParameters(4);
        Lofl = my_param.matParameters(5);
        smax = my_param.matParameters(6);
        lambda = my_param.matParameters(7);
        k = my_param.matParameters(8);
    case 'tendon material'
        g1 = my_param.matParameters(1);
        g2 = my_param.matParameters(2);
        l1 = my_param.matParameters(3);
        l2 = my_param.matParameters(4);
        lambda = my_param.matParameters(5);
        k = my_param.matParameters(6);
    case 'ogden material'
        c1=my_param.matParameters(1); %Shear-modulus-like parameter
        m1=my_param.matParameters(2); %Material parameter setting degree of non-linearity
        ksi=c1*10; %Fiber "modulus"
        alphaPar=my_param.matParameters(3);
        beta=my_param.matParameters(4);
        k_factor=my_param.matParameters(5); %Bulk modulus factor
        k=0.5.*(c1+ksi)*k_factor; %Bulk modulus
    case 'neo-Hookean fiber reinforced'
        c1=my_param.matParameters(1); %Shear-modulus-like parameter
        ksi=my_param.matParameters(2); %Fiber "modulus"
        alphaPar=my_param.matParameters(3);
        beta=my_param.matParameters(4);
        k_factor=my_param.matParameters(5); %Bulk modulus factor
        k=0.5.*(c1+ksi)*k_factor; %Bulk modulus
end

%% FEA control settings
numTimeSteps=23; %Number of time steps desired
max_refs=40; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=10; %Optimum number of iterations
max_retries=5; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=1/numTimeSteps; %Maximum time step size
symmetric_stiffness=0;

runMode=my_param.runMode; %'internal', 'external';

%% Contact parameters
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

%Contact parameters
contactInitialOffset=0.1;
contactPenalty=25;
laugon=0;
minaug=1;
maxaug=10;
fric_coeff=1e8;

%% Joining node sets
V=[V1;V2;]; %Combined node sets
E2=E2+size(V1,1); %Fixed element indices

%% Define contact surfaces

% The rigid surface of the sphere
F_contact_secondary=E2;

% The deformable surface of the slab
F_contact_primary=Fb1(Cb1==6,:);

%% Define boundary conditions

%Supported nodes
bcSupportList=unique(Fb1(Cb1==5,:)); %Bottom face
%Symmetry nodes
x_symmetryList = unique(Fb1(Cb1==3,:)); %Left face (normal to x)
y_symmetryList = unique(Fb1(Cb1==1,:)); %Front face (normal to y)

%% Define must points
timeMust = my_param.timeMust;
valMust  = dtmax.*ones(size(timeMust));
mustPoints = [timeMust valMust];

%% Defining the FEBio input structure
% See also |febioStructTemplate| and |febioStruct2xml| and the FEBio user
% manual.
%% General
%Get a template with default settings
[febio_spec]=febioStructTemplate;

%febio_spec version
febio_spec.ATTR.version='4.0';

%Module section
febio_spec.Module.ATTR.type='solid';

%% Control section
febio_spec.Control.analysis='STATIC';
febio_spec.Control.time_steps=numTimeSteps;
febio_spec.Control.step_size=1/numTimeSteps;
febio_spec.Control.solver.max_refs=max_refs;
febio_spec.Control.solver.qn_method.max_ups=max_ups;
febio_spec.Control.solver.symmetric_stiffness=symmetric_stiffness;
febio_spec.Control.time_stepper.dtmin=dtmin;
febio_spec.Control.time_stepper.dtmax=dtmax;
febio_spec.Control.time_stepper=rmfield(febio_spec.Control.time_stepper,'dtmax'); %remove default
febio_spec.Control.time_stepper.dtmax.VAL=1; 
febio_spec.Control.time_stepper.dtmax.ATTR.lc=2;
febio_spec.Control.time_stepper.max_retries=max_retries;
febio_spec.Control.time_stepper.opt_iter=opt_iter;

febio_spec.Control.output_level='OUTPUT_MUST_POINTS';

%% Material section
materialName1='Material1';
febio_spec.Material.material{1}.ATTR.name=materialName1;

switch my_param.mat_type
    case 'trans iso Mooney-Rivlin' % Transversly Isotropic Mooney-Rivlin
        k_tru=0.5*c1*k; %Bulk modulus = (initial shear modulus)x(k_factor)
        febio_spec.Material.material{1}.ATTR.type='trans iso Mooney-Rivlin';
        febio_spec.Material.material{1}.ATTR.id=1;
        febio_spec.Material.material{1}.c1=c1;
        febio_spec.Material.material{1}.c2=c2;
        febio_spec.Material.material{1}.c3=c3;
        febio_spec.Material.material{1}.c4=c4;
        febio_spec.Material.material{1}.c5=c5;
        febio_spec.Material.material{1}.lam_max=lam_max;
        febio_spec.Material.material{1}.k=k_tru;
        febio_spec.Material.material{1}.fiber.type='vector';
        febio_spec.Material.material{1}.fiber=[0,1,0];
    case 'trans iso Veronda-Westmann' % Transversly Isotropic Verdona-Westmann
        k_tru=k; %Bulk modulus = (initial shear modulus)x(k_factor)
        febio_spec.Material.material{1}.ATTR.type='trans iso Veronda-Westmann';
        febio_spec.Material.material{1}.ATTR.id=1;
        febio_spec.Material.material{1}.c1=c1;
        febio_spec.Material.material{1}.c2=c2;
        febio_spec.Material.material{1}.c3=c3;
        febio_spec.Material.material{1}.c4=c4;
        febio_spec.Material.material{1}.c5=c5;
        febio_spec.Material.material{1}.lam_max=lam_max;
        febio_spec.Material.material{1}.k=k_tru;
        febio_spec.Material.material{1}.fiber.type='vector';
        febio_spec.Material.material{1}.fiber=[0,1,0];
    case 'ogden material' % Ogden with fibers
        febio_spec.Material.material{1}.ATTR.type='solid mixture';
        febio_spec.Material.material{1}.ATTR.id=1;

        %Solid component
        febio_spec.Material.material{1}.solid{1}.ATTR.type='Ogden unconstrained';
        febio_spec.Material.material{1}.solid{1}.c1=c1;
        febio_spec.Material.material{1}.solid{1}.m1=m1;
        febio_spec.Material.material{1}.solid{1}.c2=c1;
        febio_spec.Material.material{1}.solid{1}.m2=-m1;
        febio_spec.Material.material{1}.solid{1}.cp=k;

        %The passive fiber component
        febio_spec.Material.material{1}.solid{2}.ATTR.type='fiber-exp-pow';
        febio_spec.Material.material{1}.solid{2}.ksi=ksi;
        febio_spec.Material.material{1}.solid{2}.alpha=alphaPar;
        febio_spec.Material.material{1}.solid{2}.beta=beta;
        febio_spec.Material.material{1}.solid{2}.fiber.ATTR.type='vector';
        febio_spec.Material.material{1}.solid{2}.fiber.VAL=[0 1 0];

    case 'neo-Hookean fiber reinforced' % Ogden with fibers
        febio_spec.Material.material{1}.ATTR.type='solid mixture';
        febio_spec.Material.material{1}.ATTR.id=1;

        %Solid component
        febio_spec.Material.material{1}.solid{1}.ATTR.type='Mooney-Rivlin';
        febio_spec.Material.material{1}.solid{1}.c1=c1;
        febio_spec.Material.material{1}.solid{1}.c2=0;
        febio_spec.Material.material{1}.solid{1}.k=k;

        %The passive fiber component
        febio_spec.Material.material{1}.solid{2}.ATTR.type='fiber-exp-pow';
        febio_spec.Material.material{1}.solid{2}.ksi=ksi;
        febio_spec.Material.material{1}.solid{2}.alpha=alphaPar;
        febio_spec.Material.material{1}.solid{2}.beta=beta;
        febio_spec.Material.material{1}.solid{2}.fiber.ATTR.type='vector';
        febio_spec.Material.material{1}.solid{2}.fiber.VAL=[0 1 0];
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
febio_spec.Mesh.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Mesh.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
partName1='Part1';
febio_spec.Mesh.Elements{1}.ATTR.name=partName1; %Name of this part
febio_spec.Mesh.Elements{1}.ATTR.type=elementType; %Element type
febio_spec.Mesh.Elements{1}.elem.ATTR.id=(1:1:size(E1,1))'; %Element id's
febio_spec.Mesh.Elements{1}.elem.VAL=E1; %The element matrix

partName2='Part2';
febio_spec.Mesh.Elements{2}.ATTR.name=partName2; %Name of this part
febio_spec.Mesh.Elements{2}.ATTR.type='tri3'; %Element type
febio_spec.Mesh.Elements{2}.elem.ATTR.id=size(E1,1)+(1:1:size(E2,1))'; %Element id's
febio_spec.Mesh.Elements{2}.elem.VAL=E2; %The element matrix

% -> NodeSets
nodeSetName1='bcSupportList';
febio_spec.Mesh.NodeSet{1}.ATTR.name=nodeSetName1;
febio_spec.Mesh.NodeSet{1}.VAL=mrow(bcSupportList);

nodeSetName2='x_symmetryList';
febio_spec.Mesh.NodeSet{2}.ATTR.name=nodeSetName2;
febio_spec.Mesh.NodeSet{2}.VAL=mrow(x_symmetryList(:));

nodeSetName3='y_symmetryList';
febio_spec.Mesh.NodeSet{3}.ATTR.name=nodeSetName3;
febio_spec.Mesh.NodeSet{3}.VAL=mrow(y_symmetryList(:));

%% MeshDomains section
febio_spec.MeshDomains.SolidDomain.ATTR.name=partName1;
febio_spec.MeshDomains.SolidDomain.ATTR.mat=materialName1;

febio_spec.MeshDomains.ShellDomain.ATTR.name=partName2;
febio_spec.MeshDomains.ShellDomain.ATTR.mat=materialName2;

% -> Surfaces
surfaceName1='contactSurface1';
febio_spec.Mesh.Surface{1}.ATTR.name=surfaceName1;
febio_spec.Mesh.Surface{1}.quad4.ATTR.id=(1:1:size(F_contact_primary,1))';
febio_spec.Mesh.Surface{1}.quad4.VAL=F_contact_primary;

surfaceName2='contactSurface2';
febio_spec.Mesh.Surface{2}.ATTR.name=surfaceName2;
febio_spec.Mesh.Surface{2}.tri3.ATTR.id=(1:1:size(F_contact_secondary,1))';
febio_spec.Mesh.Surface{2}.tri3.VAL=F_contact_secondary;

% -> Surface pairs
contactPairName='Contact1';
febio_spec.Mesh.SurfacePair{1}.ATTR.name=contactPairName;
febio_spec.Mesh.SurfacePair{1}.primary=surfaceName1;
febio_spec.Mesh.SurfacePair{1}.secondary=surfaceName2;

%% Boundary condition section
%-> Fix boundary conditions

febio_spec.Boundary.bc{1}.ATTR.name='zero_displacement_xyz';
febio_spec.Boundary.bc{1}.ATTR.type='zero displacement';
febio_spec.Boundary.bc{1}.ATTR.node_set=nodeSetName1;
%fix bottom face to all movement
febio_spec.Boundary.bc{1}.x_dof=1;
febio_spec.Boundary.bc{1}.y_dof=1;
febio_spec.Boundary.bc{1}.z_dof=1;

febio_spec.Boundary.bc{2}.ATTR.name='zero_displacement_xyz';
febio_spec.Boundary.bc{2}.ATTR.type='zero displacement';
febio_spec.Boundary.bc{2}.ATTR.node_set=nodeSetName2;
%fix x symmetry face to all movement along y axis
febio_spec.Boundary.bc{2}.x_dof=0;
febio_spec.Boundary.bc{2}.y_dof=1;
febio_spec.Boundary.bc{2}.z_dof=0;

febio_spec.Boundary.bc{3}.ATTR.name='zero_displacement_xyz';
febio_spec.Boundary.bc{3}.ATTR.type='zero displacement';
febio_spec.Boundary.bc{3}.ATTR.node_set=nodeSetName3;
%fix y symmetry face to all movement along x axis
febio_spec.Boundary.bc{3}.x_dof=1;
febio_spec.Boundary.bc{3}.y_dof=0;
febio_spec.Boundary.bc{3}.z_dof=0;


%% Rigid section
% -> Prescribed rigid body boundary conditions
febio_spec.Rigid.rigid_bc{1}.ATTR.name='RigidFix';
febio_spec.Rigid.rigid_bc{1}.ATTR.type='rigid_fixed';
febio_spec.Rigid.rigid_bc{1}.rb=2;
febio_spec.Rigid.rigid_bc{1}.Rx_dof=1;
febio_spec.Rigid.rigid_bc{1}.Ry_dof=1;
febio_spec.Rigid.rigid_bc{1}.Rz_dof=0;
febio_spec.Rigid.rigid_bc{1}.Ru_dof=1;
febio_spec.Rigid.rigid_bc{1}.Rv_dof=1;
febio_spec.Rigid.rigid_bc{1}.Rw_dof=1;

febio_spec.Rigid.rigid_bc{2}.ATTR.name='RigidPrescribe';
febio_spec.Rigid.rigid_bc{2}.ATTR.type='rigid_displacement';
febio_spec.Rigid.rigid_bc{2}.rb=2;
febio_spec.Rigid.rigid_bc{2}.dof='z';
febio_spec.Rigid.rigid_bc{2}.value.ATTR.lc=1;
febio_spec.Rigid.rigid_bc{2}.value.VAL=-(sphereDisplacement+contactInitialOffset);
febio_spec.Rigid.rigid_bc{2}.relative=0;

%Contact section
febio_spec.Contact.contact{1}.ATTR.type=contactType;
febio_spec.Contact.contact{1}.ATTR.surface_pair=contactPairName;
febio_spec.Contact.contact{1}.two_pass=0;
febio_spec.Contact.contact{1}.laugon=laugon;
febio_spec.Contact.contact{1}.tolerance=0.2;
febio_spec.Contact.contact{1}.gaptol=0;
febio_spec.Contact.contact{1}.minaug=minaug;
febio_spec.Contact.contact{1}.maxaug=maxaug;
febio_spec.Contact.contact{1}.search_tol=0.01;
febio_spec.Contact.contact{1}.search_radius=0.1*sqrt(sum((max(V,[],1)-min(V,[],1)).^2,2));
febio_spec.Contact.contact{1}.symmetric_stiffness=0;
febio_spec.Contact.contact{1}.auto_penalty=1;
febio_spec.Contact.contact{1}.update_penalty=1;
febio_spec.Contact.contact{1}.penalty=contactPenalty;
febio_spec.Contact.contact{1}.fric_coeff=fric_coeff;

%% LoadData section
% -> load_controller
febio_spec.LoadData.load_controller{1}.ATTR.name='LC_1';
febio_spec.LoadData.load_controller{1}.ATTR.id=1;
febio_spec.LoadData.load_controller{1}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{1}.interpolate='LINEAR';
%febio_spec.LoadData.load_controller{1}.extend='CONSTANT';
febio_spec.LoadData.load_controller{1}.points.pt.VAL=[0 0; 1 1];

febio_spec.LoadData.load_controller{2}.ATTR.name='LC_2';
febio_spec.LoadData.load_controller{2}.ATTR.id=2;
febio_spec.LoadData.load_controller{2}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{2}.interpolate='LINEAR';
%febio_spec.LoadData.load_controller{2}.extend='CONSTANT';
febio_spec.LoadData.load_controller{2}.points.pt.VAL=mustPoints;

%% Output section
% -> log file
% nodal displacements
febio_spec.Output.logfile.ATTR.file=febioLogFileName;
febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';
febio_spec.Output.logfile.node_data{1}.VAL=unique(F_contact_primary)';

% % element strain energy density
% febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_sed;
% febio_spec.Output.logfile.element_data{1}.ATTR.data='sed';
% febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';

% nodal coordinates
febio_spec.Output.logfile.node_data{2}.ATTR.file=febioLogFileName_pos;
febio_spec.Output.logfile.node_data{2}.ATTR.data='x;y;z';
febio_spec.Output.logfile.node_data{2}.ATTR.delim=',';
febio_spec.Output.logfile.node_data{2}.VAL=unique(F_contact_primary)';

% rigid body data (moments and force)
febio_spec.Output.logfile.rigid_body_data{1}.ATTR.file=febioLogFileName_rigidBody;
febio_spec.Output.logfile.rigid_body_data{1}.ATTR.data='Fx;Fy;Fz;My;Mz;Mx;z';
febio_spec.Output.logfile.rigid_body_data{1}.ATTR.delim=',';

% % Elemental stresses
% febio_spec.Output.logfile.element_data{3}.ATTR.file=febioLogFileName_stress;
% febio_spec.Output.logfile.element_data{3}.ATTR.data='s1;s2;s3';
% febio_spec.Output.logfile.element_data{3}.ATTR.delim=',';
% febio_spec.Output.logfile.element_data{3}.VAL=1:size(MeshGeometry.Specimen.elements,1);

% % Elemental strains
% febio_spec.Output.logfile.element_data{2}.ATTR.file=febioLogFileName_strain;
% febio_spec.Output.logfile.element_data{2}.ATTR.data='E1;E2;E3;Exy;J';
% febio_spec.Output.logfile.element_data{2}.ATTR.delim=',';
% febio_spec.Output.logfile.element_data{2}.VAL=1:size(MeshGeometry.Specimen.elements,1);

%% Exporting the FEBio input file
% Exporting the febio_spec structure to an FEBio input file is done using
% the |febioStruct2xml| function.

febioStruct2xml(febio_spec,febioFebFileName); %Exporting to file and domNode

%% Running the FEBio analysis
% To run the analysis defined by the created FEBio input file the
% |runMonitorFEBio| function is used. The input for this function is a
% structure defining job settings e.g. the FEBio input file name. The
% optional output runFlag informs the user if the analysis was run
% succesfully.

febioAnalysis.run_filename=febioFebFileName; %The input file name
febioAnalysis.run_logname=febioLogFileName; %The name for the log file
febioAnalysis.disp_on=disp_on; %Display information on the command window
febioAnalysis.disp_log_on=disp_on; %Display convergence information in the command window
febioAnalysis.runMode=runMode;
febioAnalysis.t_check=0.25; %Time for checking log file (dont set too small)
febioAnalysis.maxtpi=1e99; %Max analysis time
febioAnalysis.maxLogCheckTime=10; %Max log file checking time

[runFlag]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!

% if runFlag==1 %i.e. a succesful run
%     %%
%     sprintf('Test number %d successful', my_param.test_ind);
%     %     % Importing nodal displacements from a log file
%     %     dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_disp),1,1);
%     %
%     %     %Access data
%     %     N_disp_mat=dataStruct.data; %Displacement
%     %     timeVec=dataStruct.time; %Time
%     %
%     %     %Create deformed coordinate set
%     %     V_DEF=N_disp_mat+repmat(V,[1 1 size(N_disp_mat,3)]);
%     %
%     %     %%
%     %     % Importing rigid body data from a log file
%     %     dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_rigidBody)); %Nodal displacements
%     %     RigidBody_mat = dataStruct.data;
%     %     RigidBody_mat = squeeze(RigidBody_mat)';
%     %     RigidBody_mat(:,1) = []; %remove Indices
%     %     RigidBody_mat = [zeros(1,size(RigidBody_mat,2)); RigidBody_mat]; % Add zero state row
%     %     RigidBody_displacement=abs(RigidBody_mat(:,3)-center_of_mass(3));
%     %     RigidBody_force_magnitude=sqrt(sum(RigidBody_mat(:,1:3).^2,2));
%     %     RigidBody_Z_force_magnitude=abs(RigidBody_mat(:,3));
% 
% end

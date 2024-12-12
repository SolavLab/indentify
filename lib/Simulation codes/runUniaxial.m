%% NOTICE: the following code is an adaptation of DEMO_febio_0001_cube_uniaxial
%Copyright (C) 2006-2023 Kevin Mattheus Moerman and the GIBBON contributors,
%taken from the GIBBON Toolbox (www.gibboncode.org) under the license
%provided therein:
%(https://github.com/gibbonCode/GIBBON/blob/master/LICENSE).

%This function runs a compression simulation on a cylinder and saves the following data at set time steps alone:
% Compression force
% Displacement and position of all nodes on the curved surface

function [febio_spec,febioAnalysis, runFlag] = runUniaxial(my_param,disp_on)

%% Control parameters
% Path names
savePath = my_param.savePath;
if ~exist(savePath,"dir")
    mkdir(savePath) % create
end
% Defining file names
febioFebFileNamePart='tempModel';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=[febioFebFileNamePart,'.txt']; %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_pos=[febioFebFileNamePart,'_pos_out.txt']; %Log file name for exporting stress sigma_z
febioLogFileName_force=[febioFebFileNamePart,'_force_out.txt']; %Log file name for exporting force

%% Specimen parameters
%Access elements, nodes, and faces from the structure
MeshGeometry=my_param.MeshGeometry;
E=MeshGeometry.Specimen.elements; %The elements
V=MeshGeometry.Specimen.nodes; %The nodes (vertices)
Fb=MeshGeometry.Specimen.facesBoundary; %The boundary faces
Cb=MeshGeometry.Specimen.boundaryMarker; %The "colors" or labels for the boundary faces
elementType = MeshGeometry.Specimen.elementType; % hex20 / hex8

%Define applied displacement 
appliedStretch = my_param.appliedStretch;
loadingOption = my_param.loadingOption; % 'compression' or 'tension'
switch loadingOption
    case 'compression'
        displacementMagnitude=-appliedStretch; %The applied stretch for uniaxial loading
    case 'tension'
        displacementMagnitude=appliedStretch; %The applied stretch for uniaxial loading
end

%% Material parameter set
%Each material type will have a different amount of parameters to consider.
switch my_param.mat_type
    case 'neo-Hookean'
        E_youngs1 = my_param.matParameters(1);
        nu1 = my_param.matParameters(2);
    case 'Ogden 1st order'
        c1 = my_param.matParameters(1);
        m1 = my_param.matParameters(2);
        k = my_param.matParameters(3);
    case 'Ogden 2nd order'
        c1 = my_param.matParameters(1);
        m1 = my_param.matParameters(2);
        c2 = my_param.matParameters(3);
        m2 = my_param.matParameters(4);
        k = my_param.matParameters(5);
    case 'Mooney-Rivlin'
        c1 = my_param.matParameters(1);
        c2 = my_param.matParameters(2);
        k = my_param.matParameters(3);
end

%% FEA control settings
numTimeSteps=20; %Number of time steps desired
max_refs=25; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=10; %Optimum number of iterations
max_retries=5; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=1/numTimeSteps; %Maximum time step size
symmetric_stiffness=0;

runMode=my_param.runMode; %'internal', 'external';

%% Defining the boundary conditions
% The visualization of the model boundary shows colors for each side of the
% cube. These labels can be used to define boundary conditions. 

%Define surface node set
logicFace=Cb==0; %Logic for current face set
Fr=Fb(logicFace,:); %The current face set
surfaceNodeList=unique(Fr(:)); %Node set part of selected face

%Define supported node set
logicFace=Cb==1; %Logic for current face set
Fr=Fb(logicFace,:); %The current face set
bcSupportList=unique(Fr(:)); %Node set part of selected face

%Prescribed displacement nodes
logicPrescribe=Cb==2; %Logic for current face set
Fr=Fb(logicPrescribe,:); %The current face set
bcPrescribeList=unique(Fr(:)); %Node set part of selected face

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
febio_spec.Control.time_stepper.max_retries=max_retries;
febio_spec.Control.time_stepper.opt_iter=opt_iter;
febio_spec.Control.time_stepper=rmfield(febio_spec.Control.time_stepper,'dtmax'); %remove default
febio_spec.Control.time_stepper.dtmax.VAL=1; 
febio_spec.Control.time_stepper.dtmax.ATTR.lc=2;

febio_spec.Control.output_level='OUTPUT_MUST_POINTS';

%% Material section
materialName1='Material1';
febio_spec.Material.material{1}.ATTR.name=materialName1;

switch my_param.mat_type
    case 'Ogden 1st order' % Ogden 1st order
        febio_spec.Material.material{1}.ATTR.type='Ogden';
        febio_spec.Material.material{1}.ATTR.id=1;
        febio_spec.Material.material{1}.c1=c1;
        febio_spec.Material.material{1}.m1=m1;
        febio_spec.Material.material{1}.k=k;
    case 'Ogden 2nd order' % Ogden 2nd order (symmetric)
        k_tru=2*c1*k; %Bulk modulus = (initial shear modulus)x(k_factor)
        febio_spec.Material.material{1}.ATTR.type='Ogden';
        febio_spec.Material.material{1}.ATTR.id=1;
        febio_spec.Material.material{1}.c1=c1;
        febio_spec.Material.material{1}.m1=m1;
        febio_spec.Material.material{1}.c2=c2;
        febio_spec.Material.material{1}.m2=m2;
        febio_spec.Material.material{1}.k=k_tru;
    case 'Mooney-Rivlin' % Mooney Rivlin
        k_tru=2*c1*k; %Bulk modulus = (initial shear modulus)x(k_factor)
        febio_spec.Material.material{1}.ATTR.type='Mooney-Rivlin';
        febio_spec.Material.material{1}.ATTR.id=1;
        febio_spec.Material.material{1}.c1=c1;
        febio_spec.Material.material{1}.c2=c2;
        febio_spec.Material.material{1}.k=k_tru;
    case 'neo-Hookean' % Neo-Hookean (MR)
        febio_spec.Material.material{1}.ATTR.type='neo-Hookean';
        febio_spec.Material.material{1}.ATTR.id=1;
        febio_spec.Material.material{1}.E=E_youngs1;
        febio_spec.Material.material{1}.v=nu1;

end

%% Mesh section
% -> Nodes
febio_spec.Mesh.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Mesh.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
partName1='Part1';
febio_spec.Mesh.Elements{1}.ATTR.name=partName1; %Name of this part
febio_spec.Mesh.Elements{1}.ATTR.type=elementType; %Element type
febio_spec.Mesh.Elements{1}.elem.ATTR.id=(1:1:size(E,1))'; %Element id's
febio_spec.Mesh.Elements{1}.elem.VAL=E; %The element matrix


% -> NodeSets
nodeSetName1='bcSupportList';
nodeSetName2='bcPrescribeList';

febio_spec.Mesh.NodeSet{1}.ATTR.name=nodeSetName1;
febio_spec.Mesh.NodeSet{1}.VAL=mrow(bcSupportList);

febio_spec.Mesh.NodeSet{2}.ATTR.name=nodeSetName2;
febio_spec.Mesh.NodeSet{2}.VAL=mrow(bcPrescribeList);

%% MeshDomains section
febio_spec.MeshDomains.SolidDomain.ATTR.name=partName1;
febio_spec.MeshDomains.SolidDomain.ATTR.mat=materialName1;


%% Boundary condition section
% -> Fix boundary conditions
febio_spec.Boundary.bc{1}.ATTR.name='zero_displacement_z';
febio_spec.Boundary.bc{1}.ATTR.type='zero displacement';
febio_spec.Boundary.bc{1}.ATTR.node_set=nodeSetName1;
febio_spec.Boundary.bc{1}.x_dof=0;
febio_spec.Boundary.bc{1}.y_dof=0;
febio_spec.Boundary.bc{1}.z_dof=1;

febio_spec.Boundary.bc{2}.ATTR.name='prescibed_displacement_z';
febio_spec.Boundary.bc{2}.ATTR.type='prescribed displacement';
febio_spec.Boundary.bc{2}.ATTR.node_set=nodeSetName2;
febio_spec.Boundary.bc{2}.dof='z';
febio_spec.Boundary.bc{2}.value.ATTR.lc=1;
febio_spec.Boundary.bc{2}.value.VAL=displacementMagnitude;
febio_spec.Boundary.bc{2}.relative=0;

if strcmp(loadingOption,'compression')
    febio_spec.Boundary.bc{3}.ATTR.name='zero_displacement_xy';
    febio_spec.Boundary.bc{3}.ATTR.type='zero displacement';
    febio_spec.Boundary.bc{3}.ATTR.node_set=nodeSetName2;
    febio_spec.Boundary.bc{3}.x_dof=1;
    febio_spec.Boundary.bc{3}.y_dof=1;
    febio_spec.Boundary.bc{3}.z_dof=0;
end

%% LoadData section
% -> load_controller
febio_spec.LoadData.load_controller{1}.ATTR.name='LC_1';
febio_spec.LoadData.load_controller{1}.ATTR.id=1;
febio_spec.LoadData.load_controller{1}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{1}.interpolate='LINEAR';
febio_spec.LoadData.load_controller{1}.points.pt.VAL=[0 0; 1 1];

febio_spec.LoadData.load_controller{2}.ATTR.name='LC_2';
febio_spec.LoadData.load_controller{2}.ATTR.id=2;
febio_spec.LoadData.load_controller{2}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{2}.interpolate='LINEAR';
febio_spec.LoadData.load_controller{2}.points.pt.VAL=mustPoints;

%% Output section
% -> log file
febio_spec.Output.logfile.ATTR.file=febioLogFileName;
febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';
febio_spec.Output.logfile.node_data{1}.VAL=unique(surfaceNodeList)';

febio_spec.Output.logfile.node_data{2}.ATTR.file=febioLogFileName_force;
febio_spec.Output.logfile.node_data{2}.ATTR.data='Rx;Ry;Rz';
febio_spec.Output.logfile.node_data{2}.ATTR.delim=',';
febio_spec.Output.logfile.node_data{2}.VAL=unique(surfaceNodeList)';

febio_spec.Output.logfile.node_data{3}.ATTR.file=febioLogFileName_pos;
febio_spec.Output.logfile.node_data{3}.ATTR.data='x;y;z';
febio_spec.Output.logfile.node_data{3}.ATTR.delim=',';
febio_spec.Output.logfile.node_data{3}.VAL=unique(surfaceNodeList)';

% Plotfile section
febio_spec.Output.plotfile.compression=0;

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

end

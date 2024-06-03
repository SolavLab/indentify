%%

clear; close all; clc;

%%
% Plot settings
fontSize=20;
faceAlpha1=0.8;
faceAlpha2=1;
edgeColor=0.25*ones(1,3);
edgeWidth=1.5;
markerSize=25;
markerSize2=50;
lineWidth=5;
lineWidth2=3;
cMap=viridis(20);

% Define analysis settings
objectiveWeights = [0.31 0.69 0 0];
Ef = 0.1; %force measurement error (normalized)
E_disp = 0.05; % displacement measurement error (normalized)

%% LOAD EXPERIMENTAL DATA
default_running_folder = getDefaultRunPath();
exp_data_type = questdlg('Data Type','Data Type','Simulation Data','Experimental Data', 'Experimental Data');
if isempty(exp_data_type)
    error('exp_data_type was left unassigned')
end
switch exp_data_type
    case 'Simulation Data'
        %TEMP_SIMULATION DATA
        fprintf('\n Select the test_data.mat of a set of synthetic test results\n\n******************\n\n');
        [~,runPath] = uigetfile(default_running_folder,'Select test_data.mat');
        if runPath == 0
            error('runPath was left unassigned')
        end
        fprintf('******************\n Loading Data...\n******************\n');
        load(fullfile(runPath,'test_data.mat'))
        temp_ref_test = [];
        fprintf('\n Select the job representing the synthetic test results\n\n******************\n\n');
        runPath = erase(runPath,'\analysis\');
        selpath = uigetdir(runPath);
        selpath = erase(selpath,runPath); %folder name
        ref_ind = str2double(selpath(7:end));
        ref_test = test{ref_ind};

        expResults = ref_test;
        force_exp = ref_test.indenter_RB_out.Fz.data;
        depth_exp = ref_test.indenter_RB_out.z.data;
        depth_exp(2:end) = -(depth_exp(2:end)-ref_test.MeshGeometry.Indenter.center_of_mass(3));
        timeMust = depth_exp / depth_exp(end);

        n = size(ref_test.pos_out.ind, 1); % Determine the number of nodes
        defaultNodeList = true(1, n); % Default nodeList as a logical array of ones
        nodeList = defaultNodeList;

    case 'Experimental Data'
        % IMPORT DIC DATA AS RETRIEVED FROM iFEA_barycentric_coordinates
        fprintf('\n Select the expResults.mat of the test results as retrieved from iFEA_barycentric_coordinates.m\n\n******************\n\n');
        [~,runPath] = uigetfile(default_running_folder,'Select expResults.mat');
        if runPath == 0
            error('runPath was left unassigned')
        end
        load(fullfile(runPath,'expResults.mat'))
        % IMPORT FORCE DATA AS VECTOR
        force_exp = [0  0  9 24 42 63 89 116 144 176 207 240 274 311 346.5 384 420.5 457 493 532 567 606 644 680 718];
        force_exp = -force_exp*9.8 / 4; % Convert [g] to [mN], and account for a quadrant of the simulation
        % IMPORT INDENTATION DEPTH DATA AS VECTOR (helps define must times)
        depth_exp = [0 0.2 0.5  1 1.5  2 2.5  3 3.5  4 4.5  5 5.5  6 6.5  7 7.5  8 8.5  9 9.5 10 10.5 11 11.5];
        timeMust = depth_exp / depth_exp(end);
end


toleranceObjectiveValue = objectiveWeights*[Ef E_disp E_disp E_disp].^2'; %cutoff range
% Material Parameters
mat_type = 'trans iso Mooney-Rivlin'; % 'trans iso Mooney-Rivlin','trans iso Veronda-Westmann','muscle material','tendon material','ogden material'
%Initial material parameter set
matParameters.c1 = 50;
matParameters.c2 = 0;
matParameters.c3 = 0;
matParameters.c4 = 0.8;
matParameters.c5 = 35;
matParameters.lam_max = 1;
matParameters.k = 1e3;

parNamesToVary = {'c1','c5'};

%Sphere parameters
numRefineStepsSphere=2;
sphereRadius=9.53/2;

%% Control Parameters
runMode = 'external'; % FEBio run mode - 'external', 'internal'
% select analysis type (currently only indentation is implemented)
analysis_type = questdlg('Analysis type','Analysis type','Indentation','Tension', 'Compression', 'Tension');
if isempty(analysis_type)
    error('analysis_type was left unassigned')
end

% Retrieve/Assign default run path for indetify's calculations
default_running_folder = getDefaultRunPath();
% Specify runPath (directory for simulation files and subfolders)
fprintf('\n Select running folder to save temp file in\n\n******************\n\n');
runPath = uigetdir(default_running_folder,'Select Running Folder');
if runPath == 0
    error('runPath was left unassigned')
end

%% Simulation Parameters
elementType = 'hex8'; % 'hex8','hex20'
%Contact parameters
contactInitialOffset=0.01;
contactPenalty=100;
fric_coeff=0.25;
laugon=0;
minaug=1;
maxaug=10;

%% Creating model geometry and mesh

load("indentify\lib\Axisymmetric Indentation\coarse_trueSize.mat")
% Offset Box
V(:,3)=V(:,3)-max(V(:,3)); %Box Z location 

%Convert elements to faces
[F,~]=element2patch(E,[],'hex8');

%Find boundary faces
[indFree]=freeBoundaryPatch(F);
Fb=F(indFree,:);

%Create faceBoundaryMarkers based on normals
[N]=patchNormal(Fb,V); %N.B. Change of convention changes meaning of front, top etc.

faceBoundaryMarker=zeros(size(Fb,1),1);

faceBoundaryMarker(N(:,1)<-0.5)=1; %Left
faceBoundaryMarker(N(:,1)>0.5)=2; %Right
faceBoundaryMarker(N(:,2)<-0.5)=3; %Front
faceBoundaryMarker(N(:,2)>0.5)=4; %Back
faceBoundaryMarker(N(:,3)<-0.5)=5; %Bottom
faceBoundaryMarker(N(:,3)>0.5)=6; %Top

meshStruct.nodes=V;
meshStruct.facesBoundary=Fb;
meshStruct.boundaryMarker=faceBoundaryMarker;
meshStruct.faces=F;
meshStruct.elements=E;
meshStruct.elementMaterialID=ones(size(E,1),1);
meshStruct.faceMaterialID=ones(size(meshStruct.faces,1),1);
sampleHeight=max(V(:,3))-min(V(:,3));
sphereDisplacement=depth_exp(end);

MeshGeometry.Specimen = meshStruct;

%% Creating triangulated sphere surface model

[E2,V2,~]=geoSphere(numRefineStepsSphere,sphereRadius);
%Offset indentor
minZ=min(V2(:,3));
V2(:,3)=V2(:,3)-minZ+max(V(:,3))+contactInitialOffset; %Sphere Z location
center_of_mass=mean(V2,1);
MeshGeometry.Indenter.elements = E2;
MeshGeometry.Indenter.nodes = V2;
MeshGeometry.Indenter.center_of_mass=mean(V2,1);
MeshGeometry.Indenter.radius = sphereRadius;

%% Simulation setup and execution
run_log.metadata.start_time_raw = now;
run_log.metadata.start_time = datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z'));
full_time = tic;
nParameters = length(fieldnames(matParameters));
parValuesIni = zeros(1,nParameters);
par_names=fieldnames(matParameters);
for i_parameter = 1:nParameters
    parValuesIni(i_parameter) = matParameters.(par_names{i_parameter});
end

analysis.mat_type = mat_type;
analysis.matParameters = parValuesIni;
analysis.MeshGeometry = MeshGeometry;
analysis.MeshGeometry.Specimen.elementType = elementType;
analysis.timeMust = timeMust';
analysis.sphereDisplacement = sphereDisplacement;
% Current file name and save path
modelName = 'tempModel';
savePath = fullfile(runPath,modelName);
analysis.savePath = savePath;
analysis.runMode = runMode;
% Start measuring elapsed time
tic
% Send (my_param,modelName,savePath) to appropriate
% GIBBON constructor and execution function
switch analysis_type
    case 'Tension'
    case 'Compression'
    case 'Indentation'
        [febio_spec,febioAnalysis,runFlag] = runAnisotropicIndentation(analysis,1);
end

analysis.runFlag = runFlag;
[~,analysis.model_name,~] = fileparts(febioAnalysis.run_logname);
analysis = getLogfileNames(analysis,febio_spec);

if analysis.runFlag == 1
    analysis = loadDataFiles(analysis);

    % Visualize force-depth curve
    cFigure; hold on;
    title('Indenter Force curves optimisation','FontSize',fontSize);
    xlabel('Indenter Depth [%]','FontSize',fontSize); ylabel('Measured Force [N]','FontSize',fontSize); hold on;

    Hn(1)=plot(100*analysis.timeMust,abs(force_exp)/1000,'k-','lineWidth',lineWidth);
    view(2); axis tight;  grid on; axis square; axis manual;
    Hn(2)=plot(100*analysis.timeMust,abs(analysis.indenter_RB_out.Fz.data),'r.-','lineWidth',lineWidth2,'markerSize',markerSize2);
    legend(Hn,{'Experiment','Simulation'},'Location','northwest');
    set(gca,'FontSize',fontSize);
    drawnow;
end

%% Create structures for optimization

[~,parIndicesToVary] = ismember(parNamesToVary,par_names);
parValuesToVary = parValuesIni(parIndicesToVary);
% Material structure
mat_struct.par_names=par_names; %Parameter names
mat_struct.par_values=parValuesIni; %Parameter values
mat_struct.par_vary_idx = parIndicesToVary; %Parameter indices

%What should be known to the objective function:
objectiveStruct.h=Hn(2);
objectiveStruct.force_exp = force_exp;
[~,pos_data,~] = getNPosMat(analysis);
objectiveStruct.pos_data = pos_data;
objectiveStruct.disp_exp = expResults.disp_out;
objectiveStruct.nodeList = nodeList;
% objectiveStruct.strain_exp = expResults.strain;
objectiveStruct.indenterRadius = sphereRadius;
objectiveStruct.objectiveWeights = objectiveWeights;
objectiveStruct.febioAnalysis=analysis;
% objectiveStruct.febioFebFileName=febioFebFileName;
objectiveStruct.mat_struct=mat_struct;
objectiveStruct.parNormFactors=parValuesToVary; %This will normalize the parameters to ones(size(P))
objectiveStruct.Pb_struct.xx_c=parValuesToVary; %Parameter constraining centre
objectiveStruct.Pb_struct.xxlim=[parValuesToVary(1)/100 parValuesToVary(1)*10;...
    parValuesToVary(2)/100     80     ]; %Parameter bounds


%Optimisation settings
maxNumberIterations=20; %Maximum number of optimization iterations
maxNumberFunctionEvaluations=maxNumberIterations*10; %Maximum number of function evaluations, N.B. multiple evaluations are used per iteration
functionTolerance=toleranceObjectiveValue; %Tolerance on objective function value
parameterTolerance=1e-6; %Tolerance on parameter variation
displayTypeIterations='iter';

objectiveStruct.method=2;

%File names of output files
% output_names.stress=fullfile(savePath,febioLogFileName_stress);
% output_names.stretch=fullfile(savePath,febioLogFileName_stretch);
% objectiveStruct.run_output_names=output_names;

%% start optimization

Pn=parValuesToVary./objectiveStruct.parNormFactors;

switch objectiveStruct.method
    case 1 %fminsearch and Nelder-Mead
        OPT_options = optimset('fminsearch'); % 'Nelder-Mead simplex direct search'
        OPT_options = optimset(OPT_options,'MaxFunEvals',maxNumberFunctionEvaluations,...
            'MaxIter',maxNumberIterations,...
            'TolFun',functionTolerance,...
            'TolX',parameterTolerance,...
            'Display',displayTypeIterations,...
            'FinDiffRelStep',1e-2,...
            'DiffMaxChange',0.5);
        [Pn_opt,OPT_out.fval,OPT_out.exitflag,OPT_out.output]= fminsearch(@(Pn) objectiveFunctionIFEA(Pn,objectiveStruct),Pn,OPT_options);
    case 2 %lsqnonlin and Levenberg-Marquardt
        OPT_options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt');
        OPT_options = optimoptions(OPT_options,'MaxFunEvals',maxNumberFunctionEvaluations,...
            'MaxIter',maxNumberIterations,...
            'TolFun',functionTolerance,...
            'TolX',parameterTolerance,...
            'Display',displayTypeIterations,...
            'FinDiffRelStep',1e-2,...
            'DiffMaxChange',0.5);
        [Pn_opt,OPT_out.resnorm,OPT_out.residual]= lsqnonlin(@(Pn) objectiveFunctionIFEA(Pn,objectiveStruct),Pn,[],[],OPT_options);
end

%%
[Fopt,OPT_stats_out]=objectiveFunctionIFEA(Pn_opt,objectiveStruct);

%%

function [Fopt,OPT_stats_out]=objectiveFunctionIFEA(Pn,objectiveStruct)

%%

analysis = objectiveStruct.febioAnalysis;
objectiveWeights = objectiveStruct.objectiveWeights;

%% Unnormalize and constrain parameters

P=Pn.*objectiveStruct.parNormFactors; %Scale back, undo normalization
P_in=P; %Proposed P

%Constraining parameters
for q=1:1:numel(P)
    [P(q)]=boxconstrain(P(q),objectiveStruct.Pb_struct.xxlim(q,1),objectiveStruct.Pb_struct.xxlim(q,2),objectiveStruct.Pb_struct.xx_c(q));
end

%% Setting material parameters

%Acces material parameters
mat_struct=objectiveStruct.mat_struct;
parIndicesToVary = mat_struct.par_vary_idx;
parValuesNow = mat_struct.par_values;
parValuesNow(parIndicesToVary)=P;

disp('SETTING MATERIAL PARAMETERS...');
disp(['Proposed (norm.): ',sprintf(repmat('%6.16e ',[1,numel(Pn)]),Pn)]);
disp(['Proposed        : ',sprintf(repmat('%6.16e ',[1,numel(P_in)]),P_in)]);
disp(['Set (constr.)   : ',sprintf(repmat('%6.16e ',[1,numel(P)]),P)]);


analysis.matParameters = parValuesNow;

disp('Done')

%% START FEBio

[febio_spec,febioAnalysis,runFlag] = runAnisotropicIndentation(analysis,0);
analysis.runFlag = runFlag;
[~,analysis.model_name,~] = fileparts(febioAnalysis.run_logname);
analysis = getLogfileNames(analysis,febio_spec);
%pause(0.1);

if runFlag==1
    % Importing analysis data
    analysis = loadDataFiles(analysis);

    if ~isempty(objectiveStruct.h)
        objectiveStruct.h.YData=abs(analysis.indenter_RB_out.Fz.data)/1000;
        drawnow;
    end

    %Derive Fopt
    obj_fun_val = calcObjFun(analysis,objectiveStruct);
    Fforce = obj_fun_val.Ff;
    Fdisp_x = obj_fun_val.Fu_x;
    Fdisp_y = obj_fun_val.Fu_y;
    Fdisp_z = obj_fun_val.Fu_z;
    FDev = objectiveWeights(1)*Fforce+...
        objectiveWeights(2)*Fdisp_x+...
        objectiveWeights(3)*Fdisp_y+...
        objectiveWeights(4)*Fdisp_z;

    switch objectiveStruct.method
        case 1
            Fopt=sum((FDev).^2); %Sum of squared differences
        case 2
            Fopt=FDev(:);%(stressDev).^2; %Squared differences
    end

    OPT_stats_out.obj_fun_val=obj_fun_val;
    OPT_stats_out.FDev=FDev;
    OPT_stats_out.Fopt=Fopt;
    OPT_stats_out.P_opt=P;

else %Output NaN
    switch objectiveStruct.method
        case 1
            Fopt=NaN;
        case 2
            Fopt=NaN(size(FDev)); %Squared differences
    end
    OPT_stats_out=[];
end

end
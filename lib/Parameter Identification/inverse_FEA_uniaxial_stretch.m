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

% Change default axes fonts.
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize',fontSize)
set(0,'defaulttextinterpreter','latex');

% Define analysis settings
objectiveWeights = [1 0];


loadingOption = 'compression'; % 'compression', 'tension', 'joint'

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
        force_exp = sum(ref_test.force_out.Rz.data,1);
        depth_exp = min(ref_test.disp_out.uz.data,[],1);
        timeMust = ref_test.disp_out.time;

        % Specimen parameters 
        cylLength = ref_test.MeshGeometry.Specimen.cylLength; % specimen length (mm)
        cylRadius = ref_test.MeshGeometry.Specimen.cylRadius; % specimen radius (mm)
        mesh_refinement_factor = ref_test.MeshGeometry.Specimen.meshf; % Mesh refinement factor, N (scalar/vector)
        appliedStretch = ref_test.appliedStretch; % (mm)

        pointSpacing=4/mesh_refinement_factor*ones(1,2); %Desired point spacing between nodes
        [meshStruct] = hexMeshCylinder(cylRadius,cylLength,pointSpacing);
        V = meshStruct.nodes;

        n = size(V, 1); % Determine the number of nodes
        nodeList = false(1, n); % Default nodeList as a logical array of ones
        nodeList(ref_test.pos_out.ind) = true;

    case 'Experimental Data'
        % IMPORT DIC DATA AS RETRIEVED FROM iFEA_barycentric_coordinates
        fprintf('\n Select the .mat files of the test results as retrieved from iFEA_barycentric_coordinates.m\n\n******************\n\n');
        [file,runPath] = uigetfile('*.mat', 'Select Experimental Results Files', 'MultiSelect', 'on');
        if runPath == 0
            error('runPath was left unassigned')
        end
        for i = 1:length(file)
            load(fullfile(runPath,file{i}))
        end
        %CHANGE TO MATCH EXPERIMENTAL DATA
%         force_exp = -weight_data'*9.80665 / 4; % Convert [g] to [mN], and account for quarter of problem
%         depth_exp = -indenter_depth;
%         timeMust = depth_exp / depth_exp(end);

end


% toleranceObjectiveValue = objectiveWeights*[Ef E_disp E_disp E_disp].^2'; %cutoff range
toleranceObjectiveValue = 0.3*1e-2;
% Material Parameters
mat_type = 'Mooney-Rivlin'; % 'trans iso Mooney-Rivlin','trans iso Veronda-Westmann','muscle material','tendon material','ogden material'
%Initial material parameter set
matParameters.c1 = 5;
matParameters.c2 = 1;
% matParameters.c3 = 0;
% matParameters.c4 = 0;
% matParameters.c5 = 0.2;
% matParameters.lam_max = 1;
matParameters.k = 1e3;
par_names=fieldnames(matParameters);

parNamesToVary = {'c1','c2'};
[~,parIndicesToVary] = ismember(parNamesToVary,par_names);

%% Model Parameters
if strcmp(loadingOption,'joint')
    loadingOption = 'compression';
    switch mat_type
        case 'Mooney-Rivlin'
            objectiveStruct.analytical = @(x) c1*x+c_2*x.^2;
        case 'Neo-hookean'
            objectiveStruct.analytical = @(x) c1*x;
    end
end

%% Control Parameters
runMode = 'external'; % FEBio run mode - 'external', 'internal'

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

pointSpacing=4/mesh_refinement_factor*ones(1,2); %Desired point spacing between nodes
[meshStruct] = hexMeshCylinder(cylRadius,cylLength,pointSpacing);
V = meshStruct.nodes;
V(:,3) = V(:,3)-min(V(:,3)); % Move center to 0
meshStruct.nodes = V;

MeshGeometry.Specimen = meshStruct;


%% Simulation setup and execution
run_log.metadata.start_time_raw = now;
run_log.metadata.start_time = datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z'));
full_time = tic;
nParameters = length(fieldnames(matParameters));
parValuesIni = zeros(1,nParameters);

for i_parameter = 1:nParameters
    parValuesIni(i_parameter) = matParameters.(par_names{i_parameter});
end

analysis.mat_type = mat_type;
analysis.matParameters = parValuesIni;
analysis.MeshGeometry = MeshGeometry;
analysis.MeshGeometry.Specimen.elementType = elementType;
analysis.timeMust = timeMust;
analysis.appliedStretch = appliedStretch;
analysis.loadingOption = loadingOption;

% Current file name and save path
modelName = 'tempModel';
savePath = fullfile(runPath,modelName);
analysis.savePath = savePath;
analysis.runMode = runMode;
% Start measuring elapsed time
tic
% Send (my_param,modelName,savePath) to appropriate
% GIBBON constructor and execution function

[febio_spec,febioAnalysis,runFlag] = runUniaxial(analysis,1);

analysis.runFlag = runFlag;
[~,analysis.model_name,~] = fileparts(febioAnalysis.run_logname);
analysis = getLogfileNames(analysis,febio_spec);

if analysis.runFlag == 1
    analysis = loadDataFiles(analysis);
    surfaceFb = meshStruct.facesBoundary(meshStruct.boundaryMarker==0,:);
    surface_nodes = unique(surfaceFb); %only contact surface nodes
    surface_nodes_inROI = find(nodeList);
    exp_symmetry_nodes = find(abs(V(surface_nodes_inROI,1))<1e-10);
    sim_symmetry_nodes = find(abs(V(surface_nodes,1))<1e-10);
    
    % Access data
    exp_def_curve_z = expResults.pos_out.z.data(exp_symmetry_nodes,end);
    sim_def_curve_z = analysis.pos_out.z.data(sim_symmetry_nodes,end);
    exp_def_curve_y = expResults.pos_out.y.data(exp_symmetry_nodes,end);
    sim_def_curve_y = analysis.pos_out.y.data(sim_symmetry_nodes,end);
    force_sim = sum(analysis.force_out.Rz.data,1);
    
    hf = cFigure; %Open figure  
%     title('Indenter Force curves optimisation','FontSize',fontSize);
    % Visualize force-depth curve
    subplot(2,2,[1 3]); hold on;
    title('Force curves optimisation','FontSize',fontSize);
    xlabel('Displacement [mm]','FontSize',fontSize); ylabel('Measured Force [N]','FontSize',fontSize, Interpreter='latex'); hold on;
    Hf(1)=plot(depth_exp,abs(force_exp)/1000,'ko','lineWidth',lineWidth);
    view(2); axis tight;  grid on; axis square; axis manual;
    Hf(2)=plot(depth_exp,abs(force_sim)/1000,'r.-','lineWidth',lineWidth2,'markerSize',markerSize2);
    legend(Hf,{'Experiment','Simulation'},'Location','northwest');
    set(gca,'FontSize',fontSize);

    subplot(2,2,2); hold on; 
    % Visualize deformed surface
    title('Specimen surface curves optimisation','FontSize',fontSize);
    xlabel('Y [mm]','FontSize',fontSize); ylabel('Z [mm]','FontSize',fontSize); hold on;
    Hd(1)=plot(exp_def_curve_y,exp_def_curve_z,'k.','lineWidth',lineWidth,'markerSize',markerSize2);
    view(2); axis tight;  grid on; axis equal;
    Hd(2)=plot(sim_def_curve_y,sim_def_curve_z,'r.','lineWidth',lineWidth2,'markerSize',markerSize);
    legend(Hd,{'Experiment','Simulation'},'Location','southeast');
    set(gca,'FontSize',fontSize);
    
    if length(parNamesToVary) == 2
        subplot(2,2,4); hold on; 
        % Visualize convergence
        title('Objective Function Space','FontSize',fontSize);
        xlabel(parNamesToVary{1},'FontSize',fontSize); ylabel(parNamesToVary{2},'FontSize',fontSize); hold on;
        Hc=scatter(parValuesIni(parIndicesToVary(1)),parValuesIni(parIndicesToVary(2)),markerSize2,'filled');
        Hc.CData = 0;
        xlim([1.5 15]);
        ylim([0 10]);
        grid on; colorbar; clim([0 0.25]);
        set(gca,'FontSize',fontSize);
    end
    drawnow;
end

%% Create structures for optimization

parValuesToVary = parValuesIni(parIndicesToVary);
% Material structure
mat_struct.par_names=par_names; %Parameter names
mat_struct.par_values=parValuesIni; %Parameter values
mat_struct.par_vary_idx = parIndicesToVary; %Parameter indices

%What should be known to the objective function:
objectiveStruct.Hf=Hf(2); %Force plot
objectiveStruct.Hd=Hd(2); %Deformation plot
objectiveStruct.Hc=Hc; %Convergence plot

objectiveStruct.sim_symmetry_nodes = sim_symmetry_nodes;
objectiveStruct.force_exp = force_exp;
[~,pos_data,~] = getNPosMat(analysis);
objectiveStruct.pos_data = pos_data;
objectiveStruct.disp_exp = expResults.disp_out;
objectiveStruct.nodeList = nodeList(unique(surfaceFb)); %only contact surface nodes
% objectiveStruct.strain_exp = expResults.strain;
objectiveStruct.objectiveWeights = objectiveWeights;
objectiveStruct.febioAnalysis=analysis;
% objectiveStruct.febioFebFileName=febioFebFileName;
objectiveStruct.mat_struct=mat_struct;
objectiveStruct.parNormFactors=parValuesToVary; %This will normalize the parameters to ones(size(P))
objectiveStruct.Pb_struct.xx_c=parValuesToVary; %Parameter constraining centre
% objectiveStruct.Pb_struct.xxlim=[parValuesToVary(1)/100 parValuesToVary(1)*10;...
%     parValuesToVary(2)/100     80     ]; %Parameter bounds
objectiveStruct.Pb_struct.xxlim=[1.5 15;...
    0     10     ]; %Parameter bounds


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
if isequal(objectiveWeights, [1, 0])
    useForceOnly = true;
else
    useForceOnly = false;
end

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

[febio_spec,febioAnalysis,runFlag] = runUniaxial(analysis,0);
analysis.runFlag = runFlag;
[~,analysis.model_name,~] = fileparts(febioAnalysis.run_logname);
analysis = getLogfileNames(analysis,febio_spec);
%pause(0.1);
timeMust = analysis.timeMust;
FDev = zeros(1,length(timeMust));

if runFlag==1
    % Importing analysis data
    analysis = loadDataFiles(analysis);
    
    %Derive Fopt
    obj_fun_val = calcObjFun_uniaxial_compr(analysis,objectiveStruct, useForceOnly);
    Fforce = obj_fun_val.Ff;
    Fdisp_r = obj_fun_val.Fu_r;
    FDev = objectiveWeights(1)*Fforce+...
        (1-objectiveWeights(1))*Fdisp_r;
    Fopt=sum((FDev).^2);

%     if isfield(objectiveStruct,'analytical')
%         stress_strain = objectiveStruct.analytical;
%     end

    if ~isempty(objectiveStruct.Hf)
        objectiveStruct.Hf.YData=abs(sum(analysis.force_out.Rz.data,1))/1000; % Kpa*mm^2 = 0.001 N
        drawnow;
    end
    
    if ~isempty(objectiveStruct.Hd)
        objectiveStruct.Hd.YData=analysis.pos_out.z.data(objectiveStruct.sim_symmetry_nodes,end); 
        objectiveStruct.Hd.XData=analysis.pos_out.y.data(objectiveStruct.sim_symmetry_nodes,end); 
        drawnow;
    end

    if ~isempty(objectiveStruct.Hc)
        objectiveStruct.Hc.YData=[objectiveStruct.Hc.YData,P(2)]; 
        objectiveStruct.Hc.XData=[objectiveStruct.Hc.XData,P(1)];
        objectiveStruct.Hc.CData=[objectiveStruct.Hc.CData;Fopt];
        drawnow;
    end



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
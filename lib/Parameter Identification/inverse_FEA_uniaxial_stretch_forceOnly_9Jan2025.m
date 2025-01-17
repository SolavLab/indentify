%%

clear; close all; clc;

%%
% Plot settings
fontSize=20;
markerSize=25;
markerSize2=50;
lineWidth=5;
lineWidth2=3;

% Change default axes fonts.
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize',fontSize)
set(0,'defaulttextinterpreter','latex');

%%

% Define analysis settings
cylRadius = 25/2;
cylLength = 34/2; %half for symmetry
mesh_refinement_factor = 1;

initial_area = cylRadius^2*pi; %mm

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
        force_exp_comp = abs(sum(ref_test.force_out.Rz.data,1));
        if strcmp(loadingOption,'joint') || strcmp(loadingOption,'compression')
            depth_exp = min(ref_test.disp_out.uz.data,[],1);
        else
            depth_exp = max(ref_test.disp_out.uz.data,[],1);
        end
        timeMust = ref_test.disp_out.time;

        % Specimen parameters 
        cylLength = ref_test.MeshGeometry.Specimen.cylLength; % specimen length (mm)
        cylRadius = ref_test.MeshGeometry.Specimen.cylRadius; % specimen radius (mm)
        mesh_refinement_factor = ref_test.MeshGeometry.Specimen.meshf; % Mesh refinement factor, N (scalar/vector)
        appliedStretch = ref_test.appliedStretch; % (mm)

    case 'Experimental Data'
        % IMPORT DIC DATA AS RETRIEVED FROM iFEA_barycentric_coordinates
        fprintf('\n Select the .txt file of the test results of force and compression\n\n******************\n\n');
        [file,runPath] = uigetfile('*.txt', 'Select Experimental Results File');
        if runPath == 0
            error('runPath was left unassigned')
        end
        file = fullfile(runPath, file);
        % Read the entire file into a table
        data = readtable(file, 'Delimiter', '\t','ConsecutiveDelimitersRule' ,'join',...
            'VariableNamingRule','preserve');

        % Extract the weight excluding missing snapshots
        force_data = data.("Force (N)");

        % Extract the depth excluding missing snapshots
        disp_data = data.("Position (mm)")/2; %fit to half a cylinder
        %CHANGE TO MATCH EXPERIMENTAL DATA
        force_exp_comp = abs(force_data') * 1000; % Convert [N] to [mN]
        timeMust = disp_data / disp_data(end);
        appliedStretch = abs(disp_data(end));
        depth_exp = abs(disp_data);

        if strcmp(loadingOption,'joint')
            [file,runPath]=uigetfile('*.csv','Select uniaxial tension DIC data');
            if runPath == 0
                error('runPath was left unassigned')
            end
            f = fullfile(runPath, file);
            [T,p_names,vars_unique, delimiter] = loadTable(f);

            % Arrange data
            [im] = getSnapshot(T,vars_unique, delimiter); % DIC data structure ogranized according to individual images.
            [P] = getDataPoints(im,p_names,2); % Clean NaN points and reorganize

            nPoints = numel(P);
            nImages = numel(im);

            DIC3D.strain   = zeros(nImages, 1, nPoints);

            DIC3D.strain(:,1,:) = reshapeFieldData(P, 'EI', nImages, nPoints);
            DIC3D.strain(:,2,:) = reshapeFieldData(P, 'EII', nImages, nPoints);
            DIC3D.strain(:,3,:) = reshapeFieldData(P, 'strain_vonMises', nImages, nPoints);


            [file,runPath] = uigetfile('*.txt','Select results data',runPath);
            if runPath == 0
                error('runPath was left unassigned')
            end
            file = fullfile(runPath, file);
            % Read the entire file into a table
            data = readtable(file, 'Delimiter', '\t','ConsecutiveDelimitersRule' ,'join',...
                'VariableNamingRule','preserve');
            
            % Extract the weight excluding missing snapshots
            force_exp_tens = data{[im.File]+1,"Force (N)"};

            axial_strain = zeros(1,nImages);
            transverse_strain = axial_strain;

            for i=1:nImages %Retrieve shape of cylinder for each time step
                transverse_strain(i)=mean(DIC3D.strain(i,2,:),'omitnan'); %normalized transverse stretch
                axial_strain(i)=real(mean(DIC3D.strain(i,1,:),'omitnan')); %normalized axial stretch
            end

            axial_stretch = sqrt(2*axial_strain + 1);
            transverse_stretch = sqrt(2*transverse_strain + 1);

            stress_exp_tens = (force_exp_tens' ./ initial_area) .* axial_stretch * 1000; %kPa

        end

end


% toleranceObjectiveValue = objectiveWeights*[Ef E_disp E_disp E_disp].^2'; %cutoff range
toleranceObjectiveValue = 1e-4;
% Material Parameters
mat_type = 'Mooney-Rivlin'; % 'trans iso Mooney-Rivlin','trans iso Veronda-Westmann','muscle material','tendon material','ogden material'
%Initial material parameter set
matParameters.c1 = 11;
matParameters.c2 = 0;
% matParameters.c3 = 0;
% matParameters.c4 = 0;
% matParameters.c5 = 0.2;
% matParameters.lam_max = 1;
matParameters.k = 1e3;
par_names=fieldnames(matParameters);

parNamesToVary = {'c1'};
[~,parIndicesToVary] = ismember(parNamesToVary,par_names);

%% Model Parameters
if strcmp(loadingOption,'joint')
    loadingOption = 'compression';
    objectiveStruct.stress_exp_tens = stress_exp_tens;
    objectiveStruct.axial_stretch = axial_stretch;
    objectiveStruct.initial_area = initial_area;
    switch mat_type
        case 'Mooney-Rivlin'
            objectiveStruct.analytical = @(x,c1,c2) 2.*(c1+c2.*x.^(-1)).*(x.^2-x.^(-1));
        case 'Neo-hookean'
            objectiveStruct.analytical = @(x,c1) 2.*c1.*(x.^2-x.^(-1));
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

    % Access data
    force_sim = sum(analysis.force_out.Rz.data,1);
    
    hf = cFigure; %Open figure  
    % Visualize force-depth curve
    subplot(2,2,[1 3]); hold on;
    title('Force curves optimisation','FontSize',fontSize);
    xlabel('Displacement [mm]','FontSize',fontSize); ylabel('Measured Force [N]','FontSize',fontSize, Interpreter='latex'); hold on;
    Hf(1)=plot(depth_exp,abs(force_exp_comp)/1000,'ko','lineWidth',lineWidth);
    view(2); axis tight;  grid on; axis square; axis manual;
    Hf(2)=plot(depth_exp,abs(force_sim)/1000,'r.-','lineWidth',lineWidth2,'markerSize',markerSize2);
    legend(Hf,{'Experiment','Simulation'},'Location','northwest');
    set(gca,'FontSize',fontSize);
    
    if length(parNamesToVary) == 2
        subplot(2,2,[2 4]); hold on; 
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
if length(parNamesToVary) == 2
    objectiveStruct.Hc=Hc; %Convergence plot
else
    objectiveStruct.Hc=[];
end

objectiveStruct.force_exp = force_exp_comp;
objectiveStruct.febioAnalysis=analysis;
objectiveStruct.mat_struct=mat_struct;
objectiveStruct.parNormFactors=parValuesToVary; %This will normalize the parameters to ones(size(P))
objectiveStruct.Pb_struct.xx_c=parValuesToVary; %Parameter constraining centre
objectiveStruct.Pb_struct.xxlim=[1.5 30]; %Parameter bounds


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
        [Pn_opt,OPT_out.resnorm,OPT_out.residual,~,~,~,OPT_out.jacobian]= lsqnonlin(@(Pn) objectiveFunctionIFEA(Pn,objectiveStruct),Pn,[],[],OPT_options);
end

%%
[Fopt,OPT_stats_out]=objectiveFunctionIFEA(Pn_opt,objectiveStruct);


%%

function [Fopt,OPT_stats_out]=objectiveFunctionIFEA(Pn,objectiveStruct)

%%

analysis = objectiveStruct.febioAnalysis;

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
    obj_fun_val = calcObjFun_uniaxial_compr(analysis,objectiveStruct, 1);
    FDev = obj_fun_val.Ff;

    if ~isempty(objectiveStruct.Hf)
        objectiveStruct.Hf.YData=abs(sum(analysis.force_out.Rz.data,1))/1000; % Kpa*mm^2 = 0.001 N
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

    if isfield(objectiveStruct,'analytical')
        stress_exp_tens = objectiveStruct.stress_exp_tens;
        axial_stretch =objectiveStruct.axial_stretch;
        analytical = objectiveStruct.analytical;

        parValuesCell = num2cell(parValuesNow(1:(end-1)));
        stress_res = analytical(axial_stretch,parValuesCell{:})-stress_exp_tens;
        stress_res_mag = abs(stress_res);
        stress_exp_mag = max(abs(stress_exp_tens));
        squared_normalized_stress_res = (stress_res_mag ./ stress_exp_mag).^2;
        squared_normalized_stress_res(isnan(squared_normalized_stress_res))=0;

        FDev = [FDev/length(FDev),squared_normalized_stress_res/length(squared_normalized_stress_res)];
        switch objectiveStruct.method
            case 1
                Fopt=sum((FDev).^2); %Sum of squared differences
            case 2
                Fopt=FDev(:);%(stressDev).^2; %Squared differences
        end
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

function pos = reshapeFieldData(P, fieldName, nImages, nPoints)
% Extract the field values from the structure array and reshape
pos = reshape([P.(fieldName)], nImages, nPoints);
end
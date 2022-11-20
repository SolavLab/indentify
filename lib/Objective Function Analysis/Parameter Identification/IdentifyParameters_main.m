close all
clear all
clc

default_running_folder = 'C:\FEBio_runs\axisymmetry\mesh_convergence_benchmark\new\temp1\test_1';

%% Specify and load baseline parameter simulation files
% Specify baselinePath (directory for simulation files)
if ~exist(default_running_folder)
    default_running_folder = getenv('DEFAULT_RUNNING_FOLDER');
end
baselinePath = uigetdir(default_running_folder,'Select baseline running folder');
if baselinePath == 0
    error('runPath was left unassigned')
end
baseline_fname = fullfile(baselinePath,'test.mat');
baseline_spec_name = fullfile(baselinePath,'febio_spec.mat');
if ~exist(baseline_fname)
    error('test.mat wasn''t found in selected directory');
end
if ~exist(baseline_spec_name)
    error('febio_spec.mat wasn''t found in selected directory');
end
baseline = load(baseline_fname).test; % load test
if ~baseline.runFlag
    error('selected baseline test simulation did not terminate successfully');
end
baseline.febio_spec = load(baseline_spec_name).febio_spec; % load febio_spec

%% read data from logfiles
baseline = loadDataFiles(baseline);

%% Optimisation settings
    maxNumberIterations=20; %Maximum number of optimization iterations
    maxNumberFunctionEvaluations=maxNumberIterations*10; %Maximum number of function evaluations, N.B. multiple evaluations are used per iteration
    functionTolerance=1e-7; %Tolerance on objective function value
    parameterTolerance=1e-4; %Tolerance on parameter variation
    displayTypeIterations='iter';
    % Initial parameters (for optimization)
    P_exp = [baseline.p1 baseline.p2];
    P = [1.8162    0.8892].*P_exp;
    eta = 0.5;
%% Organize baseline ("experimental") data
% global force_exp indentation_depth_exp pos_exp disp_exp 
force_exp = baseline.indenter_RB_out.Fz.data;
indentation_depth_exp = baseline.indenter_RB_out.z.data;
[~,pos_exp,disp_exp] = getNPosMat(baseline);
%% Define febioAnalysis structure for testing trial parameter sets
trial_dir = fullfile(baseline.savePath,'Opt'); % create subfolder for objective function evaluation trial runs
mkdir(trial_dir);
febioFebFileName = fullfile(trial_dir,strcat(baseline.model_name,'.feb'));
febioLogFileName = fullfile(trial_dir,strcat(baseline.model_name,'.txt'));
febioAnalysis.run_filename=febioFebFileName; %The input file name
febioAnalysis.run_logname=febioLogFileName; %The name for the log file
febioAnalysis.disp_on=1; %Display information on the command window
% febioAnalysis.runMode='internal';
febioAnalysis.maxLogCheckTime=300; %Max log file checking time


%% Create structures for optimization
% Material structure
% mat_struct.id=1; %Material id
% mat_struct.mat_type = baseline.mat_type;

% docNode=set_mat_par_FEBIO(FEB_struct.run_filename,FEB_struct.run_filename,{mat_struct});

febioAnalysis.disp_on=0;
febioAnalysis.disp_log_on=0;

% What should be known to the objective function:
% obj_data.hl=hl(2); % handle for simulated force-displacement curve
% obj_data.ht2=ht(3); % handle for parameter text
% obj_data.ht3=ht(4); % handle for Fopt
% obj_data.ht4=ht(5); % handle for title

obj_data.force_exp=force_exp;
obj_data.indentation_depth_exp=indentation_depth_exp;
obj_data.pos_exp=pos_exp; 
obj_data.disp_exp = disp_exp;

obj_data.febioAnalysis=febioAnalysis;
obj_data.febio_spec=baseline.febio_spec;
obj_data.febioFebFileName=febioFebFileName;
obj_data.febioModelName = baseline.model_name;
obj_data.savePath = trial_dir;
obj_data.eta = eta;

obj_data.mat_struct.mat_type=baseline.mat_type;
obj_data.mat_struct.k_factor = baseline.k_factor;
obj_data.parNormFactors=P; %This will normalize the parameters to ones(size(P))
obj_data.P_exp = P_exp;

%         obj_data.Pb_struct.xx_c=P; %Parameter constraining centre
%         obj_data.Pb_struct.xxlim=[[P(1)/100 P(2)/100 P(3)/100]' [P(1)*100 P(2)*100 P(3)*100]']; %Parameter bounds
% obj_data.indenterVertex=E2_joined(1,1);
obj_data.iterationNumber=0;

%File names of output files
% output_names.displacement=fullfile(savePath,febioLogFileName_disp);
% output_names.force=fullfile(savePath,febioLogFileName_rigidBody);
% obj_data.run_output_names=output_names;
obj_data.indenterRadius = baseline.Indenter.indenterRadius;


%% variable_tracking_struct
global variable_tracking_struct;%<-create a global for tracking data
variable_tracking_struct=struct;
variable_tracking_struct.isfirst_flg = 1;
variable_tracking_struct.Fopt_previous = 15;

%% start optimization
obj_data.formulation = 6;
maxNumberIterations = 10;
maxNumberFunctionEvaluations = 10*maxNumberIterations;
parameterTolerance =1e-2;
rng default
Pn=P./obj_data.parNormFactors;
obj_data.method=2; % 1=> fminsearch and Nelder-Mead, 2=> lsqnonlin and Levenberg-Marquardt, 3=> fminsearch and Nelder-Mead (force+displacment),4=> lsqnonlin and Levenberg-Marquardt (force+displacment)
switch obj_data.method
    case 1 %fminsearch and Nelder-Mead
        OPT_options=optimset('fminsearch'); % 'Nelder-Mead simplex direct search'
        
        OPT_options = optimset(OPT_options,'MaxFunEvals',maxNumberFunctionEvaluations,...
                                'MaxIter',maxNumberIterations,'TolFun',functionTolerance,...
                                'TolX',parameterTolerance,'Display',displayTypeIterations,...
                                'FinDiffRelStep',1e-2,'DiffMaxChange',0.5);
                            
        [Pn_opt,OPT_out.fval,OPT_out.exitflag,OPT_out.output]= fminsearch(@(Pn) iFEARoutine(Pn,obj_data),Pn,OPT_options);
        
        case 2 %lsqnonlin and Levenberg-Marquardt
            OPT_options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt');
            OPT_options = optimoptions(OPT_options,'MaxFunEvals',maxNumberFunctionEvaluations,...
                'MaxIter',maxNumberIterations,...
                'TolFun',functionTolerance,...
                'TolX',parameterTolerance,...
                'Display',displayTypeIterations,...
                'FinDiffRelStep',1e-2,...
                'DiffMaxChange',0.5);
            [Pn_opt,OPT_out.resnorm,OPT_out.residual]= lsqnonlin(@(Pn) iFEARoutine(Pn,obj_data),Pn,[],[],OPT_options);
end

%%
test_arr = variable_tracking_struct.test_arr;
for ii=1:numel(test_arr)
    p1_data(ii) = test_arr{ii}.p1;
    p2_data(ii) = test_arr{ii}.p2;
    Fopt(ii) = test_arr{ii}.obj_fun_val.Fopt;
end
plot3(p1_data,p2_data,Fopt,'-o')

%%
%  c1_data = [];
%         m1_data = [];
%         F_opt_data = [];
%         figure;
%         hold on
%         for qt=1:1:num_iterations %Loop over time increments
%             c1_data(qt) = variable_tracking_struct.P_all{qt}(1);
%             m1_data(qt) = variable_tracking_struct.P_all{qt}(2);
%             F_opt_data(qt) = variable_tracking_struct.Fopt_all{qt}(1);
%         end
     figure;
        hold on
    num_iterations=numel(test_arr);
    c1_data = p1_data;
    m1_data = p2_data
    F_opt_data = Fopt;
    c1_tru = obj_data.P_exp(1)
    m1_tru = obj_data.P_exp(2)
    
        [~, opt_ind] = min(F_opt_data);
        C_map = colormap(jet(num_iterations));
        for i=1:num_iterations-1
            plot3(c1_data(i:i+1)/c1_tru,m1_data(i:i+1)/m1_tru, F_opt_data(i:i+1), '-', 'Color', C_map(i,:), 'HandleVisibility', 'Off');
        end
        xlabel('c1_i');
        ylabel('m1_i');
        zlabel('F_{opt}_i');
        % P_opt=Pn_opt.*obj_data.parNormFactors;
        % colormap(C_map);
        scatter_size = 36;
        scatter3(c1_data/c1_tru,m1_data/m1_tru, F_opt_data,scatter_size*ones(1,num_iterations),C_map, 'HandleVisibility' , 'Off'); % All points
        scatter3(c1_data(opt_ind)/c1_tru,m1_data(opt_ind)/m1_tru, F_opt_data(opt_ind),1.5*scatter_size,C_map(end,:),'filled', 'DisplayName', 'Optimal'); %Optimal point
        scatter3(c1_data(1)/c1_tru,m1_data(1)/m1_tru, F_opt_data(1),1.5*scatter_size,C_map(1,:),'filled', 'DisplayName', 'Innitial'); %Initial point
        scatter3(c1_tru/c1_tru,m1_tru/m1_tru,0, 'x', 'DisplayName', 'Tru Values'); % tru point
        grid on;
        legend show;
        caxis([0 num_iterations])
        cbar = colorbar;
        cbar.Title.String = 'Itteration #';
        title(OPT_out.output.message);
        disp(['Optimum of objective Fcn: ', num2str(F_opt_data(opt_ind))]);
        disp(['Optimal C1: ', num2str(c1_data(opt_ind))]);
        disp(['Optimal m1: ', num2str(m1_data(opt_ind))]);
%         disp(['Optimal k: ', num2str(variable_tracking_struct.P_all{opt_ind}(3))]);
        ax1 = gca;
        ax1.ZScale = 'log';
        zlim([0, 1.01]);
        t =xline(c1_tru/c1_tru, 'k'); yline(m1_tru/m1_tru, 'k'); 

        % disp(['C1_error', num2str(abs((c1_data(opt_ind)-c1_tru)/c1_tru))]);
        % disp(['m1_error', num2str(abs((c1_data(opt_ind)-c1_tru)/c1_tru))]);
        % disp(['m1_error', num2str(abs((c1_data(opt_ind)-c1_tru)/c1_tru))]);
        figure;
        semilogy(F_opt_data);
        xlabel('k (itteration)');
        ylabel('F_{opt}^k')
%%


function [Fopt]=iFEARoutine(Pn,obj_data)
    global variable_tracking_struct
    persistent n_calls;% this will be incrementing each time the function is called
    if isempty(n_calls)
        n_calls=0;
    end
    %%
    febioFebFileName=obj_data.febioFebFileName;
    febio_spec=obj_data.febio_spec;
    savePath = fileparts(obj_data.febioFebFileName);
    febioLogFileName_rigidBody = febio_spec.Output.logfile.rigid_body_data{1}.ATTR.file;
    febioLogFileName_disp = febio_spec.Output.logfile.node_data{1}.ATTR.file;
    febioLogFileName_pos = febio_spec.Output.logfile.node_data{2}.ATTR.file;
   
    %% Unnormalize and constrain parameters
    P_trial=Pn.*obj_data.parNormFactors; %Scale back, undo normalization
    P_exp = obj_data.P_exp;
    %% print to screen
    disp('SETTING MATERIAL PARAMETERS...');
    disp(['Proposed (norm.): ',sprintf(repmat('%6.16e ',[1,numel(Pn)]),Pn)]);
    disp(['Proposed        : ',sprintf(repmat('%6.16e ',[1,numel(P_trial)]),P_trial)]);
    disp(['True   : ',sprintf(repmat('%6.16e ',[1,numel(P_exp)]),P_exp)]);
    %% Update material values in febio_spec and write to new FEB file
    mat_struct = obj_data.mat_struct;
    febio_spec.Material.material{1} = UpdateMat(febio_spec.Material.material{1},mat_struct,P_trial);
    febioStruct2xml(febio_spec,febioFebFileName); %Exporting to file and domNode
    %% setup trial_test structure
    trial_test.p1 = P_trial(1);
    trial_test.p2 = P_trial(2);
    trial_test.model_name = obj_data.febioModelName;
    trial_test.savePath = obj_data.savePath;
    %% START FEBio
    tic
        [trial_test.runFlag]=runMonitorFEBio(obj_data.febioAnalysis); %run
    trial_test.elapsed_time = toc;
    % load log data
    trial_test = getLogfileNames(trial_test,febio_spec);
    trial_test = loadDataFiles(trial_test);
    %% Evaluate Objective function
    [trial_test.obj_fun_val] = objFun(trial_test,obj_data);
    eta = obj_data.eta;
    trial_test.obj_fun_val.Ffu = simpleInterp(trial_test.obj_fun_val.Ff,trial_test.obj_fun_val.Fu,eta);
    trial_test.obj_fun_val.Ffu_avg = mean(trial_test.obj_fun_val.Ffu);
    Fopt = mean(trial_test.obj_fun_val.Ffu_avg);
    %% Store data and exit function
    % increment each time the objective_function completes
    n_calls=n_calls+1;
    % save plot variables for future animation
    variable_tracking_struct.test_arr{n_calls} = trial_test;
end


function [mat] =  UpdateMat(mat,mat_struct,p)
    mat_type = mat_struct.mat_type;
    k_factor = mat_struct.k_factor;
    switch mat_type
        case 'Ogden_1st'
            mat.c1 = p(1);
            mat.m1 = p(2);
            mu = p(1)/2; % initial shear modulus
            mat.k = k_factor*mu;
        case 'Ogden_symmetric' 
            mat.c1 = p(1);
            mat.c2 = p(2);
            mu = p(1); % initial shear modulus
            mat.k = k_factor*mu;
        case 'Mooney-Rivlin'
            mat.c1 = p(1);
            mat.c2 = p(2);
            mu = 2*(p(1)+p(2));
            mat.k = k_factor*mu; % initial shear modulus
        case  'Neo-Hookean(Young)'
            mat.E = p(1);
            mat.v = p(2);
    end
end
% Post-proccesing script for creating the objective functions contour plots
% and hessian metrics heatmaps (Figs. 4-7 in paper)
clear; close all; clc;

%% Post-processing controls  <<<<< User-specified

colormap_data_field ='Fval_contour';%'Fval_surface'; %'Fval_contour','grad_mag_contour','scatter3Fval'

% Construct cell arrays for "by exp_params_carry (must be in same
% structure)" if same analysis is to be run many times.
exp_params_carray = {}; mat_type_carray = {};
exp_params_carray{end+1}=[12.000000000000002, 4.5, 1000.0]; mat_type_carray{end+1} = 'Mooney-Rivlin';

%% Specify analysis folders, baseline parameters and material model names  <<<<< User-specified

% % % % % % % % % % % EDIT AS NEEDED % % % % % % % % % % % <<<<< User-specified (START)

% Retrieve/Assign default run path for indetify's calculations
default_running_folder = getDefaultRunPath();
% Specify runPath (directory for simulation files and subfolders)
fprintf('******************\n Specify runPath (directory for simulation files and subfolders)\n******************\n');
runPath = uigetdir(default_running_folder,'Select Running Folder');
if runPath == 0
    error('runPath was left unassigned')
end

dir_name_carray = dir(runPath);
dir_name_carray = dir_name_carray([dir_name_carray.isdir]); %remove files
dir_name_carray(ismember({dir_name_carray.name},{'.','..'})) = []; % remove dot directories
dir_name_carray = {dir_name_carray.name}'; %turn from structure into cell array
if strncmp(dir_name_carray{1},'test',4)||strcmp(dir_name_carray{1},'analysis') %if this is a one material subfolder
    dir_name_carray = {runPath}; %set run path as the directory name
else
    dir_name_carray = append(runPath,'\',dir_name_carray(:,1));
end

file_name_carray = fullfile(dir_name_carray,'test_data.mat');

%% Specify Reference (synthetic experimental data)
specify_ref_test = questdlg('Choose reference data selection method','reference data selection','by exp_params_carry (must be in same structure)','manually', 'by exp results','by exp_params_carry (must be in same structure)');

%% Set controls (overwritings)
% revaluate objective function (if exists already)
override_obj_fun_val = questdlg('Override previous objective function evaluations?','Override evaluations','yes','no', 'no');
% save objective function relvauations
if strcmp(override_obj_fun_val,'yes')
    save_obj_fun_val = questdlg('Save n objfunction evaluations?','Override evaluations','yes','no', 'no');
else
    save_obj_fun_val = 'no';
end

%% Run over each folder

for dir_ind = 1:numel(dir_name_carray)
    % Update wait bar with each passing folder
    save_obj_fun_val_temp = save_obj_fun_val;
    max_obj_fun_val = 0;
    local_runPath = dir_name_carray{dir_ind};
    % Create folder to store analysis
    dir_analysis = fullfile(local_runPath,'\analysis');
    load(fullfile(local_runPath,'run_log.mat'));
    mat_type = run_log.metadata.mat_type;
    if ~exist(dir_analysis,"dir")
        mkdir(dir_analysis);
        % Load and adjust test data. This adds the output data (FEBio logfile data) from
        % each job and unites all test data to a single structure 'test' which is
        % saved in the sub folder labeled 'analysis' as 'test_data.mat'
        test = RunLog2Test(run_log); %this adds aditional field to test (warning: may remove data);
        file_name = fullfile(dir_analysis,'test_data.mat');
        save(file_name,'test','-v7.3');
        fprintf('******************\n%s\nsuccessfully saved\n******************\n',file_name);
    else
        fprintf('******************\n Loading Data...\n******************\n');
        load(fullfile(dir_analysis,'test_data.mat'))
    end

    %% Objective Function Evaluation
    % Meta data of simulations
    nSteps = length(test{1}.timeMust)-1; % Check amount of time steps in the simulation
    X = run_log.metadata.X;

    % Create objective function structure data:
    % * For all simulations skipped for efficiency the value will be -inf by default
    % * For failed simulations the value will be NaN (will later be interpolated)
    % * For converged simulations the value will be calculated
    F = zeros(size(X))-inf; TF = zeros(size(X));
    for t=1:nSteps
        timeStep = strcat('t_',num2str(t));  % Define field name
        objectiveValues.Ff.(timeStep)=F;
        objectiveValues.Fu_r.(timeStep)=F;
    end

    % Get Experimental data
    switch specify_ref_test
        case 'by exp_params_carry (must be in same structure)' % find by exp_params_carry
            exp_params = exp_params_carray{dir_ind};
            ref_ind = [];
            % find reference test id.
            for i=1:numel(test)
                if test{i}.runFlag==0
                    continue;
                end
                if norm([test{i}.matParameters]-exp_params)<1e-10
                    ref_ind = i; %define the index of the test used as reference
                    break;
                end
            end
            if isempty(ref_ind)
                error('Reference test index not found!')
            end
            ref_test = test{ref_ind};
        case 'manually' % manually select the job representing the synthetic test results (specify subfolder of job).
            %  meshes must be identical but material models and parameters may vary.
            temp_ref_test = [];
            fprintf('\n Select the job representing the synthetic test results\n\n******************\n\n');
            selpath = uigetdir(local_runPath);
            selpath = erase(selpath,local_runPath); %folder name
            ref_ind = str2double(selpath(7:end));
            ref_test = test{ref_ind};
        case 'by exp results'
            fprintf('******************\n Select exp_results file \n******************\n');
            % WRITE CODE THAT MATCHES RESULTS TO SIMULATION DATA
    end

    % Select partial data
    cylLength = ref_test.MeshGeometry.Specimen.cylLength;
    cylRadius = ref_test.MeshGeometry.Specimen.cylRadius;
    [~,pos_data,~] = getNPosMat(ref_test);
    ROI_bounds_normalized = 0.8; %boundaries in units of cylLength
    ROI_bounds = ROI_bounds_normalized*cylLength; %absolute boundaries
    nodeList = (pos_data(:,1) <= ROI_bounds(1));
    objectiveStruct.nodeList = nodeList;

    ref_test.disp_out.ux.data = ref_test.disp_out.ux.data(nodeList,:);
    ref_test.disp_out.uy.data = ref_test.disp_out.uy.data(nodeList,:);
    ref_test.disp_out.uz.data = ref_test.disp_out.uz.data(nodeList,:);

    % Insert values to objectiveStruct (used for evaluating the objective function).
    objectiveStruct.force_exp = sum(ref_test.force_out.Rz.data,1);
    objectiveStruct.pos_data = pos_data(:,:,1);
    objectiveStruct.disp_exp = ref_test.disp_out;

    %     % Add synthetic random noise to the synthetic test data measurements
    %     S_pos = rand(size(objectiveStruct.pos_exp));
    %     S_force = rand(size(objectiveStruct.force_exp));
    %     de = 0; % de (%) <<< set de=0 to disable synthetic noise
    %     M_pos = de*((2*S_pos-1)/100)+1; % =100+-de%
    %     M_force = de*((2*S_force-1)/100)+1; % =100+-de%
    %     objectiveStruct.pos_exp = objectiveStruct.pos_exp .*M_pos;
    %     objectiveStruct.force_exp = objectiveStruct.force_exp .*M_force;


    %Calculate objective function value for each test
    for i=1:numel(test)     % run over tests
        if isfield(test{i},'obj_fun_val')
            if strcmp(override_obj_fun_val,'yes') %re-evaluate objFun
                if test{i}.runFlag==1
                    warning('Evaluating test #%d/%d in %s.',i,numel(test),dir_analysis);
                    test{i}.obj_fun_val = calcObjFun_uniaxial_compr(test{i},objectiveStruct,0); %evaluate objFun
                end
            end
        else
            warning('Obective function field not found. Evaluating test #%d/%d in %s.',i,numel(test),dir_analysis);
            if test{i}.runFlag==0 %simulation failed
                % Construct a missing value for each objective function result
                test{i}.obj_fun_val.Ff=[0,NaN(1,nSteps)];
                test{i}.obj_fun_val.Fu_r=[0,NaN(1,nSteps)];

            elseif test{i}.runFlag==2 %simulation skipped
                continue; % Skip calculations on tests that were not simulated
            else %simulation should be evaluated
                test{i}.obj_fun_val = calcObjFun_uniaxial_compr(test{i},objectiveStruct,0); %evaluate objFun
                save_obj_fun_val_temp = 'yes';
            end
        end

        % Export values for Hessian calculations
        for t=1:nSteps % Loop over the time steps
            timeStep = strcat('t_',num2str(t));  % Define field name
            objectiveValues.Ff.(timeStep)(i) = test{i}.obj_fun_val.Ff(t+1);
            objectiveValues.Fu_r.(timeStep)(i) = test{i}.obj_fun_val.Fu_r(t+1);
        end
    end
    % Check if there are any NaN values in the data, and if so interpolate
    % based on adjacent evaluated data. If the failed data is a range of
    % adjacent data, try interpolating using a different dimension.
    TF_temp = TF;
    for dim=1:ndims(X)
        F = objectiveValues.Ff.(timeStep);
        idx = find(isnan(F));
        first = true;
        if ~isempty(idx)
            for t=1:nSteps
                timeStep = strcat('t_',num2str(t));
                [objectiveValues.Ff.(timeStep),TF_temp] = fillmissing(objectiveValues.Ff.(timeStep),'linear',dim,'EndValues','nearest');
                objectiveValues.Fu_r.(timeStep) = fillmissing(objectiveValues.Fu_r.(timeStep),'linear',dim,'EndValues','nearest');
                if first
                    TF=TF_temp; % save the initial map of failed simulations
                    first=false;
                end
            end
            for i=1:length(idx) %check all interpolated data and reset
                if isinf(objectiveValues.Ff.(timeStep)(idx(i)))
                    for t=1:nSteps
                        timeStep = strcat('t_',num2str(t));
                        objectiveValues.Ff.(timeStep)(idx(i))=NaN;
                        objectiveValues.Fu_r.(timeStep)(idx(i))=NaN;
                    end
                end
            end
        end
    end

    switch save_obj_fun_val_temp
        case 'yes'
            fprintf('******************\n Saving test data...\n******************\n');
            warning('Saving test data in %s.',dir_analysis);
            save(fullfile(dir_analysis,'test_data.mat'),'test','-v7.3');
            save_obj_fun_val_temp = 'no'; %reset
    end

    %% Post Processing
    % Create finalized identifiability results

    % Determine best accuracey for Hessian calculations
    sizeX = size(X);
    n_accuracy = min(sizeX(sizeX>1))-1;
    if mod(n_accuracy,2); n_accuracy=n_accuracy-1; end %if n is odd, take the next smallest even number
    if n_accuracy>8; n_accuracy=8; end
%     n_accuracy = 2; % uncomment this line to override automatic accuracy

    % normalize parameter space
    num_dim = length(sizeX(sizeX>1)); % the number of dimensions
    % find index of Hessian center point in the parameter space
    S = cell (1,num_dim); % create a cell array to store the output arguments
    [S{:}] = ind2sub(size(X),ref_ind); % convert the linear index of the reference point to subscripts
    mid_value = X{S{:}}; % get the value of the parameter space
    mid_value(mid_value==0)=0.1; % replace any zero values with 0.1 to avoid division by zero
    % create a normalized space
    X_norm = cell(size(X));
    for i=1:numel(X)
        X_norm{i} = X{i}./mid_value;
    end

    objectiveValues.Ff.sumOfSquares=F*0;
    objectiveValues.Fu_r.sumOfSquares=F*0;

    for t=1:nSteps
        timeStep = strcat('t_',num2str(t));
        objectiveValues.Ff.sumOfSquares = objectiveValues.Ff.sumOfSquares + objectiveValues.Ff.(timeStep);
        objectiveValues.Fu_r.sumOfSquares = objectiveValues.Fu_r.sumOfSquares + objectiveValues.Fu_r.(timeStep);
    end

    % Determine objective function shapes of each test
    Hf = getHessian(X_norm,objectiveValues.Ff.sumOfSquares,n_accuracy);
    Hr = getHessian(X_norm,objectiveValues.Fu_r.sumOfSquares,n_accuracy);
    Ef = 0.05; %force measurement error (normalized)
    E_disp = 0.05; % displacement measurement error (normalized)
    p=zeros(num_dim,length(fieldnames(objectiveValues))); %initialize best objective function data matrix
    err=zeros(num_dim,1);
    std=err;
    clc;
    for i=1:num_dim
        % Define the objective function as a function handle
        s = struct ('type','()','subs',{{i,i}});
        objfun = @(x) sqrt(x*([Ef E_disp].^2'))* ...
            sqrt(2*subsref(inv(x(1)*Hf + x(2)*Hr),s));
        % Define the constraint function as a function handle
        constrfun = @(x) deal(-subsref(inv(x(1)*Hf + x(2)*Hr),s), x(1) + x(2) - 1);
        % Define an initial guess for the weights
        x0 = [0.5 0.5];
        % Define some options for the solver
        options = optimoptions('fmincon','Display','iter');
        % Define the lower and upper bounds for the weights
        lb = [0 0];
        ub = [1 1];
        % Call the fmincon function with the bounds
        [p(i,:),fval] = fmincon(objfun,x0,[],[],[],[],lb,ub,constrfun,options);
        std(i) = sqrt(p(i,:)*[Ef E_disp].^2');
%         err(i)=std(i)*sqrt(2*fval);
        err(i)=fval;
    end
    % Display the optimal results for each parameter
    disp(mat_type)
    for i=1:num_dim
        % objective function shape to use
        OF = p(i,1)*objectiveValues.Ff.sumOfSquares+p(i,2)*objectiveValues.Fu_r.sumOfSquares;
        H = getHessian(X_norm,OF,n_accuracy);
        [V,K] = eig(H);
        varNames={'Hessian','Eigenvectors','Eigenvalues'};
        timeStep = table(H,V,K,'VariableNames',varNames);
        disp(timeStep)
        fprintf('Error of "%s" is ±%d (±%.2f%%)\n',run_log.metadata.varried_parameters{i},err(i),100*err(i));
        fprintf('Where the optimal objective function uses %.2f*Ff+%.2f*Fu_r\nAnd the standard deviation is S=%.2f\n\n',p(i,:),std(i));
    end


end

%% Figure Generation <<<<< User-specified
% This code plots a 2D slice of the objective function space as a 3D
% surface and a 2D contour. The contour levels are chose to highlight the
% the confidence interval.
fontSize = 15;
% Initialize variables and parameters
close all
operateOn = [1,2]; % choose the dimensions that the local hessian will operate on!
idealParam = 1; % choose objective function shape to use
idealP = p(idealParam,:);
OF = idealP(1)*objectiveValues.Ff.sumOfSquares+idealP(2)*objectiveValues.Fu_r.sumOfSquares;
% OF = 0.5*objectiveValues.Ff.sumOfSquares+0.5*objectiveValues.Fu_r.sumOfSquares;
% get the names of the parameters
varried_parameters = run_log.metadata.varried_parameters;
fields = run_log.metadata.fields;

%determine which dimensions in space have depth
n = length(fields);
var = zeros(n,1);
var(ismember(fields, varried_parameters)) = 1;
idx1=find(var); idx0 = find(~var); %get the indices of the positions of the varried parameters and the fixed parameters
num_dim = ndims(X);
v = 1:num_dim;

% find index of Hessian center point in the parameter space
S = cell (1,num_dim);
[S{:}] = ind2sub(size(X),ref_ind); % convert the linear index of the reference point to subscripts
S=cell2mat(S);
target_ind = ones(n,1);
target_ind(var==1) = S; % update the target indices with the subscripts of the reference point for the varried parameters


% Reorder the arrays and indices according to the dimensions to operate on
% such that the first two values correspond to operateOn
temp_target = target_ind; temp_target(idx1(operateOn))=[];
temp_target = [target_ind(idx1(operateOn));temp_target];
var_temp = var; var(idx1(operateOn))=[];
var_temp = [var_temp(idx1(operateOn));var];
idx_temp=find(var_temp); % find the indices of the varried parameters in the reordered vector
local_X = permute(X_norm, [operateOn, setdiff(v, operateOn)]); %set the array such that the first and second dimensions define a matrix to operate on
local_Z = permute(OF, [operateOn, setdiff(v, operateOn)]); % set the array such that the first and second dimensions define the function values to plot
local_TF = permute(TF, [operateOn, setdiff(v, operateOn)]); % set the array such that the first and second dimensions define the logical values of the interpolated data


% Get the 2D slice of the arrays and the coordinates
C=num2cell(temp_target(idx_temp)); % convert the target indices to a cell array
C{1}=':'; C{2}=':'; %structure the input local_X(:,:,target_ind(other dimensions))
local_X = local_X(C{:}); % get 2D slice of the normalized positions
local_Z = local_Z(C{:}); % get 2D slice of the function values
idx = local_TF(C{:}); % get 2D slice of the logical values
xData=zeros(size(local_X));
yData=zeros(size(local_X));
for i=1:numel(xData) % loop over all the elements
    xData(i)=local_X{i}(idx1(operateOn(1))); % assign the x-coordinate as the value of the varried parameter corresponding to the first dimension
    yData(i)=local_X{i}(idx1(operateOn(2))); % assign the y-coordinate as the value of the varried parameter corresponding to the second dimension
end


% Plot the 2D slice as a 3D surface and a 2D contour
figure() % create a new figure
surf(xData,yData,local_Z) %draw 3d surface
shading interp % use interpolated shading
xlabel(varried_parameters{operateOn(1)}) % label the x-axis with the name of the varried parameter
ylabel(varried_parameters{operateOn(2)}) % label the y-axis with the name of the varried parameter
figure() % create another figure

% Define countour lines
flatFlag = 0;
stdDeviation = std(idealParam)^2;
% stdDeviation = 0.033^2;
numCurves = 10;
[Fx,Fy] = gradient(local_Z);
grad_temp = sqrt(Fx.^2+Fy.^2);

min_grad = min(abs(grad_temp),[],'all');
contour_lines = linspace(min(local_Z,[],'all'), max(local_Z,[],'all'), numCurves);
[~, std_idx] = min(abs(contour_lines - stdDeviation)); % Find the index of the value closest to S
[~, cont_idx] = min(abs(local_Z - 2*stdDeviation),[],'all'); % Find the index of the point closest to S
contour_lines(std_idx) = stdDeviation;
contour_lines(std_idx:end) = stdDeviation:(2*abs(grad_temp(cont_idx))):(stdDeviation+(numCurves-std_idx)*2*abs(grad_temp(cont_idx)));
contour_lines(1:std_idx-1) = linspace(min(local_Z,[],'all'),stdDeviation,std_idx-1);
if std_idx==1 contour_lines=[0,contour_lines]; std_idx=2; numCurves=numCurves+1; end
contourf(xData,yData,local_Z,contour_lines) %draw 2d contours with data points
%Have all figures with the same colorbar:
% global_max = max([objectiveValues.Ff.t_5,objectiveValues.Fu_x.t_5,objectiveValues.Fu_y.t_5,objectiveValues.Fu_z.t_5],[],'all');
% clim([0, 0.5*global_max]) % uncomment this line to set the color limits
%Optimize colorbar:
clim([0, (stdDeviation+(numCurves-std_idx)*2*abs(grad_temp(cont_idx)))])
colorbar % add a colorbar

hold on 
axis manual

% Scatter simulation run values on plot
scatter (xData,yData,'k',"+") % scatter the data points in black
if any(idx,'all') %if any point was interpolated, mark it red
    scatter (xData(idx),yData(idx), 'r','filled') % scatter the interpolated points in red
end
xlabel(append('$',varried_parameters{operateOn(1)},'/',varried_parameters{operateOn(1)},'^*$'),'FontSize',fontSize,'Interpreter','latex') % label the x-axis with the name of the varried parameter
ylabel(append('$',varried_parameters{operateOn(2)},'/',varried_parameters{operateOn(2)},'^*$'),'FontSize',fontSize,'Interpreter','latex') % label the y-axis with the name of the varried parameter

% Draw estimated ellipse
theta_fun = @(x,y) ([x;y]-[1;1]);
Hessian_temp = getHessian(local_X,local_Z,n_accuracy);
func = @(x,y) dot((theta_fun(x,y)'*Hessian_temp)',theta_fun(x,y));
fcontour(func,'--w', 'LevelList', std(idealParam)^2, 'Visible', 'on','tag','Contours');

hold off 

%% Optimal Objective Function

% % Define the objective function as a function handle
% objfun = @(x) -det(x(1)*Hf + x(2)*Hx + x(3)*Hy + x(4)*Hz);
% % Define the constraint function as a function handle
% constrfun = @(x) deal([], x(1) + x(2) + x(3) + x(4) - 1);
% % Define an initial guess for the weights
% x0 = [0.25 0.25 0.25 0.25];
% % Define some options for the solver
% options = optimoptions('fmincon','Display','iter');
% % Define the lower and upper bounds for the weights
% lb = [0 0 0 0];
% ub = [1 1 1 1];
% % Call the fmincon function with the bounds
% [p,fval] = fmincon(objfun,x0,[],[],[],[],lb,ub,constrfun,options);

% % Define the objective function as a function handle
% s = struct ('type','()','subs',{{1,1}});
% objfun = @(x) abs(subsref(inv(x(1)*Hf + x(2)*Hx + x(3)*Hy + x(4)*Hz),s));
% % Define the constraint function as a function handle
% constrfun = @(x) deal(-10*subsref(inv(x(1)*Hf + x(2)*Hx + x(3)*Hy + x(4)*Hz),s), x(1) + x(2) + x(3) + x(4) - 1);
% % Define an initial guess for the weights
% x0 = [0.25 0.25 0.25 0.25];
% % Define some options for the solver
% options = optimoptions('fmincon','Display','iter');
% % Define the lower and upper bounds for the weights
% lb = [0 0 0 0];
% ub = [1 1 1 1];
% % Call the fmincon function with the bounds
% [p,fval] = fmincon(objfun,x0,[],[],[],[],lb,ub,constrfun,options);

% Axisymmetric Indentation main script
close all
clear all
clc
%% USER-DEFINED SETTINGS
    % Material Parameters
    mat_type = 'OG'; % 'OG','OM','MR','NH'
    P1 = [1:5:51]*1e-3; % Range of first material parameter (scalar/vector)
    P2 = [1:3:37]; % Range of second material parameter (scalar/vector)
    k_factor = 1000; % Range of bulk material parameter multiplier (scalar/vector)
    % Specimen parameters
    R_sp = 60; % specimen radius (mm)
    H_sp = 60; % specimen height (mm)
    mesh_refinement_factor = 2; % Mesh refinement factor, N (scalar/vector)
    ALPHA = 2; % angle of sector
    R_bias = 0.5; % bias factor in radial direction (R_bias=(beta_r-1) from paper)
    H_bias = 0.5; % bias factor in axial direction (H_bias=(beta_r-1) from paper)
    elementType = 'hex20'; % 'hex8','hex20'
    benchmark_flag = 0; % 0 - axisymmetric model, 1 - full 3D model.
    ignore_formula = 1; % 0 - special formula for element size bias, 1 - same as in paper
    % Indenter Parameters
    R_ind = 15; % indenter radius (mm)
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
runPath = uigetdir(default_running_folder,'Select Running Folder');
if runPath == 0
    error('runPath was left unassigned')
end
%% Simulation setup and execution
run_log.metadata.start_time_raw = now;
run_log.metadata.start_time = datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z'));
full_time = tic;
test_ind = 1;
for meshf_ind = 1:length(mesh_refinement_factor)
    for P1_ind =1:length(P1)
        for P2_ind = 1:length(P2) % parfor (P2_ind = 1:length(P2),5)
            for k_factor_ind=1:length(k_factor)
                test_ind = k_factor_ind+(P2_ind-1)*length(k_factor)+(P1_ind-1)*length(P2)+(meshf_ind-1)*length(P1);
                test = [];
            % Current Material Parameters
                test.mat_type = mat_type;
                test.p1 = P1(P1_ind);
                test.p2 = P2(P2_ind);
                test.k_factor = k_factor(k_factor_ind);
            % Additional Control parameters
                test.test_ind = test_ind;
                test.P1_ind = P1_ind;
                test.P2_ind = P2_ind;
                test.k_factor_ind = k_factor_ind;
                test.Specimen.meshf = mesh_refinement_factor(meshf_ind);
                test.Specimen.alpha = ALPHA;
                test.Specimen.elementType = elementType;
                test.Specimen.benchmark_flag = benchmark_flag;
                test.Specimen.ignore_formula = ignore_formula;
                test.Specimen.NR = round(10*test.Specimen.meshf); %number of elements along radial direction
                test.Specimen.NH = round(5*test.Specimen.meshf); %number of elements along axial direction
                test.Specimen.R_bias = R_bias;
                test.Specimen.H_bias = H_bias;
                test.Specimen.R = R_sp;
                test.Specimen.H = H_sp;
                test.Indenter.indenterRadius = R_ind;
            % Progress data 
                T = table(test_ind, mesh_refinement_factor(meshf_ind));
                T.Properties.VariableNames = {'test ID', 'meshf'};
                disp(T);
            % Current file name and save path
                modelName = strcat('test_',num2str(test_ind)); %regular mesh
                savePath = fullfile(runPath,modelName);
                test.savePath = savePath;
            % Parameters of current job 
                my_param = [];
                my_param.mat_type = test.mat_type;
                my_param.p1 = test.p1;
                my_param.p2 = test.p2;
                my_param.runMode = runMode;
                my_param.k_factor = test.k_factor;
                my_param.Specimen = test.Specimen;
                my_param.Indenter = test.Indenter;
            % Start measuring elapsed time
                tic
            % Send (my_param,modelName,savePath) to appropriate
            % GIBBON constructor and execution function
                switch analysis_type
                    case 'Tension'
                    case 'Compression'
                    case 'Indentation'
                            [febio_spec,febioAnalysis, runFlag,savePath,MeshGeometry] = runIndentation(my_param,modelName,savePath);
                end
                [~,test.model_name,~] = fileparts(febioAnalysis.run_logname);
                test.node_data_files = {};
                test.element_data_files = {};
                test.rigid_body_data = {};
                test.Specimen.alpha = MeshGeometry.Specimen.alpha; %might be updated after run (adaptive)
                test.MeshGeometry = MeshGeometry;
                test = getLogfileNames(test,febio_spec);
                test.elapsed_time = toc;
                test.runFlag = runFlag;
                if ~runFlag %i.e. an unsuccesful run
                    warning(['Error termination (test ID=', num2str(test_ind), 'model name =',modelName]);
                end
                % Store test metadata in run_log struct
                if isfield(test,'MeshGeometry')
                    run_log.test_light{test_ind} = rmfield(test,'MeshGeometry');
                else
                    run_log.test_light{test_ind} = test;
                end
                run_log.test{test_ind} = test;
                parsave(fullfile(savePath,'test.mat'),'test'); % save test.mat
                % write to text (use copy withouth MeshGeometry field to
                % avoid clutter)
                if isfield(test,'MeshGeometry')
                    yaml.WriteYaml(fullfile(savePath,'test.txt'),rmfield(test,'MeshGeometry'));
                else
                    yaml.WriteYaml(fullfile(savePath,'test.txt'),test);
                end
                % save febio_spec.mat (used for parameter identification analysis)
                parsave(fullfile(savePath,'febio_spec.mat'),'febio_spec');
            end
        end
    end
end
% Update and save run_log structure
run_log.metadata.mat_type = mat_type;
run_log.metadata.P1 = P1;
run_log.metadata.P2 = P2;
run_log.metadata.mesh_refinement_factor = mesh_refinement_factor;
run_log.metadata.runPath = runPath;
run_log.metadata.end_time_raw = now;
run_log.metadata.end_time = datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z'));
run_log.metadata.run_time = sprintf('%d hours, %d minutes and %f seconds',floor(toc(full_time)/3600), rem(floor(toc(full_time)/60),60), rem(toc(full_time),60));
save(fullfile(runPath,'run_log.mat'),'run_log');
yaml.WriteYaml(fullfile(runPath,'run_log.txt'),rmfield(run_log,'test'))
% find and report failed jobs 
bad_ind = [];
for i=1:numel(run_log.test) 
    if ~run_log.test{i}.runFlag
        bad_ind(end+1) = i;
    end
end
disp(['Successful runs: ', num2str(numel(run_log.test)-length(bad_ind)), '/',num2str(numel(run_log.test))]);
disp(['Failed runs: ', num2str(length(bad_ind)), '/',num2str(numel(run_log.test))]);
% open run path in explorer
winopen(runPath);
% Load and adjust test data. This adds the output data (FEBio logfile data) from
% each job and unites all test data to a single structure 'test' which is saved in runPath
% as 'test.mat'
%
% load_arrange_flag = questdlg('Automatically load and arragne data structure?','Load and Arragne Data Structure','yes','no', 'yes');
% if strcmp(load_arrange_flag, 'yes')
    % Load run_log file (metadata of all jobs)
    test = RunLog2Test(run_log); %this adds aditional field to test (warning: may remove data);
    file_name = fullfile(runPath,'test.mat');
    save(file_name,'test','-v7.3');
    fprintf('******************\n%s\nsuccessfully saved\n******************\n',file_name);
% end
toc(full_time) % print elapsed time to command window
%% 
% _*indentify footer text*_ 
% 
% License: <https://github.com/SolavLab/indentify/blob/main/LICENSE>
% 
% indentify: An open-source project for exploring the identifiability of 
% soft-tissue material parameters from noninvasive indentation test and
% inverse finite-element analysis.
% 
% Copyright (C) 2022 Zohar Oddes, Dana Solav, and the indentify contributors
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Affero General Public License for more details.
% 
% You should have received a copy of the GNU Affero General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
clc
close all

% This function creates a run_log.mat file in the specified runPath directory, which contains the metadata and test results of each simulation.
% The function takes the same parameters that were used to create the failed file as input arguments, and then rebuilds the run_log.mat file by reading
% the YAML files from the subfolders of each simulation.
% The function also saves the multidimensional array X, which contains the parameter values for each analysis,
% in the run_log.mat file. The function returns the run_log structure as output.

%Sphere parameters
numRefineStepsSphere=2;
sphereRadius=4;
contactInitialOffset=0.01;
sampleHeight=sphereRadius*6; %Height

% Retrieve/Assign default run path for indetify's calculations
default_running_folder = getDefaultRunPath();
% Specify runPath (directory for simulation files and subfolders)
runPath = uigetdir(default_running_folder,'Select Running Folder');
if runPath == 0
    error('runPath was left unassigned')
end
dir_name = dir(runPath);
dir_name = dir_name(strncmp({dir_name.name},'test',4)); %remove files
dir_name = {dir_name.name};
dir_size = size(dir_name,2);

yaml_file = append(runPath,'\run_log.txt');
run_log = yaml.ReadYaml(yaml_file);
meta=run_log.metadata;
mat_type = meta.mat_type; % 'trans iso Mooney-Rivlin','trans iso Veronda-Westmann','muscle material','tendon material','ogden material'
paramValuesForAnalyses = meta.paramValuesForAnalyses;
paramValuesForAnalyses = cell2mat(paramValuesForAnalyses);
nAnalyses = length(paramValuesForAnalyses);
nParameters = length(meta.varried_parameters);
varried_parameters = cell(1,nParameters);
for i=1:nParameters
    varried_parameters(i) = meta.varried_parameters{i};
end
C = zeros(1,nParameters);
for i=1:size(paramValuesForAnalyses,2)
    C(i) = length(unique(paramValuesForAnalyses(:,i)));
end

X = cell(C);

%Recreate Sphere data for indentation test
[E2,V2,~]=geoSphere(numRefineStepsSphere,sphereRadius);
%Offset indentor
minZ=min(V2(:,3));
V2(:,3)=V2(:,3)-minZ+(sampleHeight/2)+contactInitialOffset; %Sphere Z location
center_of_mass=mean(V2,1);
MeshGeometry.Indenter.center_of_mass=mean(V2,1);
MeshGeometry.Indenter.radius = sphereRadius;


%Build a full run_log file of only simulations that were "skipped" as basis
total_run_time=tic;
waitbar_build = waitbar(0,' ','Name','FEBio Simulations', 'HandleVisibility', 'On');
for i=1:nAnalyses
    test.mat_type = mat_type;
    test.matParameters = paramValuesForAnalyses(i,:);
    test.runFlag = 2;
    test.test_ind = i;
    test.MeshGeometry=MeshGeometry;
    modelName = strcat('test_',num2str(i)); %regular mesh
    savePath = fullfile(runPath,modelName);
    test.savePath = savePath;
    X{i} = paramValuesForAnalyses(i,:);
    run_log.test{i}=test;
    estimate_time = seconds((nAnalyses-i)*toc(total_run_time)/i);
    estimate_time.Format = 'dd:hh:mm:ss';
    waitbar(i/nAnalyses,waitbar_build,sprintf('Filling runlog (%d of %d)\n Time remaining: About %s',i,nAnalyses,floor(estimate_time)), 'HandleVisibility', 'On')
end
%Read all YAML files and insert them into the full run_log
total_run_time=tic;
for i=1:dir_size
    yaml_file = append(runPath,'\',dir_name{i},'\test.txt');
    test = yaml.ReadYaml(yaml_file);
    test.MeshGeometry=MeshGeometry;
    run_log.test{test.test_ind} = test;
    estimate_time = seconds((dir_size-i)*toc(total_run_time)/i);
    estimate_time.Format = 'dd:hh:mm:ss';
    waitbar(i/dir_size,waitbar_build,sprintf('Reading YAML (%d of %d)\n Time remaining: About %s',i,dir_size,floor(estimate_time)), 'HandleVisibility', 'On')
end
%Build general metadata:
run_log.metadata.mat_type = mat_type;
run_log.metadata.runPath = runPath;
run_log.metadata.fields = meta.fields;
run_log.metadata.varried_parameters = varried_parameters;
run_log.metadata.X = X;
waitbar(1,waitbar_build,'Saving runlog', 'HandleVisibility', 'On')
save(fullfile(runPath,'run_log.mat'),'run_log','-v7.3');
close(waitbar_build);

% Bulk simulation main script
clear; close all; clc;

%% USER-DEFINED SETTINGS

% Material Parameters
mat_type = 'Mooney-Rivlin'; % 'trans iso Mooney-Rivlin','trans iso Veronda-Westmann','muscle material','tendon material','ogden material', 'neo-Hookean fiber reinforced'
leftSide = linspace(1, 6.535, 15 + 1); rightSide = linspace(6.535, 21, 15 + 1);
matParameters.c1 = [leftSide(1:end-1), rightSide]; % Range of first material parameter (scalar/vector)
matParameters.c2 = linspace(0,0,1); % Range of second material parameter (scalar/vector)
% matParameters.c3 = linspace(0,0,1); % Range of third material parameter (scalar/vector)
% matParameters.c4 = linspace(0,0,1); % Range of fourth material parameter (scalar/vector)
% matParameters.c5 = linspace(330.26*0.25,330.26*1.75,25); % Range of fifth material parameter (scalar/vector)
% matParameters.P6 = linspace(1,1.1,1);
% matParameters.lam_max = 1;
matParameters.k = 1e3; % Range of bulk material parameter multiplier (scalar/vector)

% Specimen geometry
mesh_path = 'C:\Users\user\OneDrive - Technion\Amit-Dana shared folder\Research\Parameter Identification\Mesh Anlysis\100x100x60 cube\mesh_2.mat';
%Indenter parameters (ignore if doing tension/compression)
numRefineStepsSphere=2;
sphereRadius=9.53/2;
sphereDisplacement=14;

% Specimen parameters (ignore if doing indentation)
cylLength = 20; % specimen length (mm)
cylRadius = 13; % specimen radius (mm)
mesh_refinement_factor = 1.5; % Mesh refinement factor, N (scalar/vector)
appliedStretch = 8; % (mm)

elementType = 'hex8'; % 'hex8','hex20'

%% Control Parameters
runMode = 'external'; % FEBio run mode - 'external', 'internal'
% select analysis type (currently only indentation is implemented)
testType = questdlg('Analysis type','Analysis type','Indentation','Tension', 'Compression', 'Tension');
if isempty(testType)
    error('analysis_type was left unassigned')
end
% decide if to run all simulations in space or just those needed to
% calcualte the Hessian at the center point of space
run_all = questdlg('Run every simulation in the parameter space?','Simulation Size','Yes','No','Yes');

% Retrieve/Assign default run path for indetify's calculations
default_running_folder = getDefaultRunPath();
% Specify runPath (directory for simulation files and subfolders)
runPath = uigetdir(default_running_folder,'Select Running Folder');
if runPath == 0
    error('runPath was left unassigned')
end

%Contact parameters
contactInitialOffset=0.01;
contactPenalty=100;
fric_coeff=0.25;
laugon=0;
minaug=1;
maxaug=10;

%% Define must points

timeMust = [0 0.2 0.4 0.6 0.8 1]'; %List of time points in which a simulation result will be available.

%% Creating model geometry and mesh
switch testType
    case 'Indentation'
        load(mesh_path)
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
        MeshGeometry.Specimen = meshStruct;

        % Creating triangulated sphere surface model
        [E2,V2,~]=geoSphere(numRefineStepsSphere,sphereRadius);
        %Offset indentor
        minZ=min(V2(:,3));
        V2(:,3)=V2(:,3)-minZ+max(V(:,3))+contactInitialOffset; %Sphere Z location
        center_of_mass=mean(V2,1);
        MeshGeometry.Indenter.elements = E2;
        MeshGeometry.Indenter.nodes = V2;
        MeshGeometry.Indenter.center_of_mass=mean(V2,1);
        MeshGeometry.Indenter.radius = sphereRadius;

    otherwise
        pointSpacing=4/mesh_refinement_factor*ones(1,2); %Desired point spacing between nodes
        [meshStruct] = hexMeshCylinder(cylRadius,cylLength,pointSpacing);
        V = meshStruct.nodes;
        V(:,3) = V(:,3)-min(V(:,3)); % Move center to 0 
        meshStruct.nodes = V;

        MeshGeometry.Specimen = meshStruct;
end



%% Simulation setup and execution
run_log.metadata.start_time_raw = now;
run_log.metadata.start_time = datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z'));
full_time = tic;
nAnalyses = 1;
% Compute number of unique analyses
fields=fieldnames(matParameters);
nParameters = length(fieldnames(matParameters));
paramValuesForAnalyses = cell(1,nParameters);
paramCount = zeros(1,nParameters);
for i_parameter = 1:nParameters
    paramValuesForAnalyses{i_parameter} = matParameters.(fields{i_parameter});
    nAnalyses = nAnalyses*numel(paramValuesForAnalyses{i_parameter});
    paramCount(i_parameter) = numel(paramValuesForAnalyses{i_parameter});
end

% Compute parameter values for each analyses
[paramValuesForAnalyses{:}] = ndgrid(paramValuesForAnalyses{:});
n = nParameters;
paramValuesForAnalyses = reshape(cat(n+1, paramValuesForAnalyses{:}),[],n); % Each row is the set of parameter of one the analyses
% Allocate memory fo the output structure
analyses = cell(nAnalyses, 1); % output structure
multi_value_param=paramCount-1;
multi_value_param=find(multi_value_param);
paramCount(paramCount==1)=[];
if length(paramCount) == 1 % fix numbering in case of edge case 
    paramCount = [paramCount, 1];
end
X = cell(paramCount);

%Determine which simulations will be run or not
switch run_all
    case 'Yes'
        run_all = true;
    case 'No'
        run_all = false;
        if length(multi_value_param)<=2
            run_all = true;
        end
end

waitbar_sim = waitbar(0,'Running simulation','Name','FEBio Simulations', 'HandleVisibility', 'On');
total_run_time=tic;
for i_test = 1:nAnalyses
    % Current Material Parameters
    test.mat_type = mat_type;
    test.matParameters = paramValuesForAnalyses(i_test,:);
    X{i_test} = paramValuesForAnalyses(i_test,:); %Build material space (for Hessian computation)
    % Additional Control parameters
    test.test_ind = i_test;
    test.MeshGeometry = MeshGeometry;
    test.MeshGeometry.Specimen.elementType = elementType;
    % Current file name and save path
    modelName = strcat('test_',num2str(i_test)); %regular mesh
    savePath = fullfile(runPath,modelName);
    test.savePath = savePath;
    test.runMode = runMode;
    test.timeMust = timeMust;
    test.sphereDisplacement = sphereDisplacement;
    test.appliedStretch = appliedStretch;
    %Check if simulation will be used in Hessian calculation for material
    %properties at excat center of parameter space
    mid = paramValuesForAnalyses(ceil(nAnalyses/2),:);
    chk = test.matParameters(multi_value_param) - mid(multi_value_param);
    use_simulation = ~all(chk) || run_all; % Return true if the parameters are different or run_all is 1
    % Start measuring elapsed time
    tic
    % Send (my_param,modelName,savePath) to appropriate
    % GIBBON constructor and execution function
    switch testType
        case 'Tension'
            test.loadingOption='tension';
            if use_simulation
                [febio_spec,febioAnalysis,runFlag] = runUniaxial(test,1);
            else
                test.runFlag = 2;
                run_log.test{i_test} = test;
            continue;
            end
            test.MeshGeometry = rmfield(test.MeshGeometry,'Specimen');
            test.MeshGeometry.Specimen.cylRadius = cylRadius;
            test.MeshGeometry.Specimen.cylLength = cylLength;
            test.MeshGeometry.Specimen.meshf = mesh_refinement_factor;
        case 'Compression'
            test.loadingOption='compression';
            if use_simulation
                [febio_spec,febioAnalysis,runFlag] = runUniaxial(test,1);
            else
                test.runFlag = 2;
                run_log.test{i_test} = test;
            continue;
            end
            test.MeshGeometry = rmfield(test.MeshGeometry,'Specimen');
            test.MeshGeometry.Specimen.cylRadius = cylRadius;
            test.MeshGeometry.Specimen.cylLength = cylLength;
            test.MeshGeometry.Specimen.meshf = mesh_refinement_factor;
        case 'Indentation'
            if use_simulation
                [febio_spec,febioAnalysis,runFlag] = runAnisotropicIndentation(test,1);
            else
                test.runFlag = 2;
                run_log.test{i_test} = test;
            continue;
            end
            test.MeshGeometry.Indenter = rmfield(test.MeshGeometry.Indenter,{'nodes';'elements'});
            test.MeshGeometry = rmfield(test.MeshGeometry,'Specimen');
    end

    [~,test.model_name,~] = fileparts(febioAnalysis.run_logname);
    test.node_data_files = {};
    test.element_data_files = {};
    test.rigid_body_data = {};
    test = getLogfileNames(test,febio_spec);
    test.elapsed_time = toc;
    test.runFlag = runFlag;
    if runFlag==1 %i.e. a succesful run
        sprintf('Test number %d successful', i_test);
    end
    if ~runFlag %i.e. an unsuccesful run
        warning(['Error termination (test ID=', num2str(i_test), 'model name =',modelName]);
    end
    run_log.test{i_test} = test;
    % write to text (use copy withouth MeshGeometry field to
    % avoid clutter)
    if isfield(test,'MeshGeometry')
        yaml.WriteYaml(fullfile(savePath,'test.txt'),rmfield(test,'MeshGeometry'));
    else
        yaml.WriteYaml(fullfile(savePath,'test.txt'),test);
    end
    % save febio_spec.mat (used for parameter identification analysis)
    parsave(fullfile(savePath,'febio_spec.mat'),'febio_spec');
    % Progress data
    estimate_time = seconds((nAnalyses-i_test)*toc(total_run_time)/i_test);
    estimate_time.Format = 'dd:hh:mm:ss';
    waitbar(i_test/nAnalyses,waitbar_sim,sprintf('Running simulation (%d of %d)\n Time remaining: About %s',i_test,nAnalyses,floor(estimate_time)), 'HandleVisibility', 'On');
    %T = table(i_test, mesh_refinement_factor);
    %T.Properties.VariableNames = {'test ID', 'meshf'};
    % disp(T);
end
waitbar(1,waitbar_sim,'Compiling data...');

%% Update and save run_log structure
% find and report failed jobs
bad_ind = [];
for i=1:numel(run_log.test)
    if ~run_log.test{i}.runFlag
        bad_ind(end+1) = i;
    end
end

run_log.metadata.mat_type = mat_type;
run_log.metadata.paramValuesForAnalyses = paramValuesForAnalyses;
% run_log.metadata.mesh_refinement_factor = mesh_refinement_factor;
run_log.metadata.runPath = runPath;
run_log.metadata.timeMust = timeMust';
run_log.metadata.end_time_raw = now;
run_log.metadata.end_time = datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z'));
run_log.metadata.run_time = sprintf('%d hours, %d minutes and %f seconds',floor(toc(full_time)/3600), rem(floor(toc(full_time)/60),60), rem(toc(full_time),60));
run_log.metadata.fields = fields;
run_log.metadata.varried_parameters = fields(multi_value_param);
run_log.metadata.failed_runs = bad_ind;

yaml.WriteYaml(fullfile(runPath,'run_log.txt'),run_log.metadata);

run_log.metadata.X = X; %add multidimensional array after .yaml
waitbar(1,waitbar_sim,'Saving data...');
save(fullfile(runPath,'run_log.mat'),'run_log');

pause(2)
disp(['Successful runs: ', num2str(numel(run_log.test)-length(bad_ind)), '/',num2str(numel(run_log.test))]);
disp(['Failed runs: ', num2str(length(bad_ind)), '/',num2str(numel(run_log.test))]);
% open run path in explorer
winopen(runPath);
toc(full_time) % print elapsed time to command window
close(waitbar_sim);
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
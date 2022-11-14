% This function reads FEBio log files and adds the data as fields to the
% 'test' matlab structure. 
% Example - 
% suppose - 
%       test{:}.node_data_files = {'tempModel_disp_out.txt'  'tempModel_pos_out.txt'}
% Then  - 
%       test = loadDataFile(test)
%  test{:} will have the fields test{:}.disp_out and test{:}.pos_out
%%%%%%%%%%%%%
% e.g. 
    % test{1}.pos_out
    % 
    % ans = 
    % 
    %   struct with fields:
    % 
    %        N: 101
    %     time: [6×1 double]
    %      ind: [101×1 double]
    %        x: [1×1 struct]
    %        y: [1×1 struct]
    %        z: [1×1 struct]"
%%%%%%%%%%%%%%%%%

% If the principle strain fields is available within the log files, 
% then this function also calculates the principle streches - L1-L3, and
% the first and second invariants of C - I1 and I2.

function [test] = loadDataFiles(varargin)
    switch numel(varargin)
        case 1
            temp_test = varargin{1};
            user_output_flg = 0; % default
        case 2
            temp_test = varargin{1};
            user_output_flg = varargin{2};
    end
    % loading test{:} with data from febio log files
    field_names = {};
    if ~iscell(temp_test) % convert to cell array if not already
        test{1} = temp_test;
    else
        test = temp_test;
    end
    if user_output_flg
        f1 = waitbar(0,'Please wait...');
    end
    for i=1:numel(test)
        if user_output_flg
            waitbar(i/numel(test),f1, ['Reading tests (',num2str(i),'/', num2str(numel(test)),')'])
        end
        if test{i}.runFlag
            % gather data from node_data_files
            if isfield(test{i},'node_data_files')
                for ii=1:numel(test{i}.node_data_files) % run on each data file
                    temp_name = test{i}.node_data_files{ii}; % temporary file name
                    [~,temp_field_name,~] = fileparts(temp_name); % temporary field name (includes model name as prefix)
                    temp_field_name(1:length(test{i}.model_name)+1) = []; % remove model name from temp_field_name
                    [test{i}.(temp_field_name)] = log2struct(fullfile(test{i}.savePath,temp_name),1,0);
                    if ~any(strcmp(field_names,temp_field_name))
                        field_names{end+1} = temp_field_name;
                    end
                end
            end
            % gather data from element_data_files
            if isfield(test{i},'element_data_files')
                for ii=1:numel(test{i}.element_data_files) 
                    temp_name = test{i}.element_data_files{ii};
                    [~,temp_field_name,~] = fileparts(temp_name);
                    temp_field_name(1:length(test{i}.model_name)+1) = [];
                    [test{i}.(temp_field_name)] = log2struct(fullfile(test{i}.savePath,temp_name),1,0);
                    if ~any(strcmp(field_names,temp_field_name))
                        field_names{end+1} = temp_field_name;
                    end
                end
            end
             % gather data from rigid_body_data_files
            if isfield(test{i},'rigid_body_data_files')
                for ii=1:numel(test{i}.rigid_body_data_files)
                    temp_name = test{i}.rigid_body_data_files{ii};
                    [~,temp_field_name,~] = fileparts(temp_name);
                    temp_field_name(1:length(test{i}.model_name)+1) = [];
                    [test{i}.(temp_field_name)] = log2struct(fullfile(test{i}.savePath,temp_name),1,0);
                    if ~any(strcmp(field_names,temp_field_name))
                        field_names{end+1} = temp_field_name;
                    end
                end
            end
            if isfield(test{i},'strain_out')
                if isfield(test{i}.strain_out,'E1')
                    test{i}.strain_out.L1 = sqrt(2*test{i}.strain_out.E1.data+1);
                end
                if isfield(test{i}.strain_out,'E2')
                    test{i}.strain_out.L2 = sqrt(2*test{i}.strain_out.E2.data+1);
                end
                if isfield(test{i}.strain_out,'E3')
                    test{i}.strain_out.L3 = sqrt(2*test{i}.strain_out.E3.data+1);
                end
                if (isfield(test{i}.strain_out,'L1'))&&(isfield(test{i}.strain_out,'L2'))&&(isfield(test{i}.strain_out,'L3'))
                    L1sq = test{i}.strain_out.L1.^2;
                    L2sq = test{i}.strain_out.L2.^2;
                    L3sq = test{i}.strain_out.L3.^2;
                    test{i}.strain_out.I1 = L1sq+L2sq+L3sq;
                    test{i}.strain_out.I2 = L1sq.*L2sq+L1sq.*L3sq+L2sq.*L3sq;
                end
            end
            if numel(test)==1
                test = test{1};
            end
                    
            if isfield(test, 'MeshGeometry')
                if isfield(test,'indenter_RB_out')
                    test.indenter_RB_out.z.data = test.indenter_RB_out.z.data-test.MeshGeometry.Indenter.center_of_mass(3);
                    test.indenter_RB_out.z.data(1) = 0;
                end
            end
                    
                    
%             % stress output 
%             [test{i}.stress_out] = log2struct(fullfile(test{i}.savePath,'tempModel_stress_out.txt'),1,0);
%             [test{i}.strain_out] = log2struct(fullfile(test{i}.savePath,'tempModel_strain_out.txt'),1,0);
%             [test{i}.disp_out] = log2struct(fullfile(test{i}.savePath,'tempModel_disp_out.txt'),1,0);
%             if isfile(fullfile(test{i}.savePath,'rigid.txt'))
%                 [test{i}.rigid_out] = log2struct(fullfile(test{i}.savePath,'rigid.txt'),1,0);
%             end
        else 
            warning(['test{', num2str(test{i}.test_ind),'}.runFlag==0. No data was loaded'])  
        end
    end
    if user_output_flg
        close(f1)
        fprintf('The following fields were added to test{i} struct:\n');
        fprintf('\t# %s\n',field_names{:});
    end
    
%     % Caulculate principal stretches from principal strains
%     f2 = waitbar(0,'Please wait...');
%     for i=1:numel(test)
%         waitbar(i/numel(test),f2, ['Caulculate principal stretches (',num2str(i),'/', num2str(numel(test)),')'])
%         if test{i}.runFlag
%             test{i}.strain_out.Lambda1 = sqrt(2*test{i}.strain_out.E1.data+1);
%             test{i}.strain_out.Lambda2 = sqrt(2*test{i}.strain_out.E2.data+1);
%             test{i}.strain_out.Lambda3 = sqrt(2*test{i}.strain_out.E3.data+1);
%          else 
%             warning(['test{', num2str(test{i}.test_ind),'}.runFlag==0. No principal streches were calculated'])  
%         end
%     end
%     close(f2)
end

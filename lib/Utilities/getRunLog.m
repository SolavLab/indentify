%when varargin is empty, the value of the enviorment variable
%DEFAULT_RUNNING_FOLDER is used
function [run_log] = getRunLog(varargin)
    if isempty(varargin)
        def_path = getenv('DEFAULT_RUNNING_FOLDER');
    else
        def_path = varargin{1};
        
    end
    runPath = uigetdir(def_path,'Select Folder');
    fname = fullfile(runPath,'run_log.mat');
    loaded_data = load(fname,'run_log');
    disp(['Successfully loaded', fname]);
    run_log = loaded_data.run_log;
    
    
    N = numel(run_log.test);
    disp(['Total number of simulations found: ', num2str(N)]); 
    disp('Batch metadata:');
    disp(struct2table(run_log.metadata));
end


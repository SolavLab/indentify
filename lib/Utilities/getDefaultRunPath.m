% This function retrieves the DEFAULT_RUNNING_FOLDER enviroment variable
% which specifies identify's default folder for carrying calculations.
% If the enviroment variable has not been defined yet then the function will
% attempt to set it according to user input. Should this attempt also fail,
% the function returnt an empty string.
function [default_running_folder] = getDefaultRunPath()
    default_running_folder = getenv('DEFAULT_RUNNING_FOLDER'); % default run path
    if isempty(default_running_folder) % check default run path
        % Set default run path
        set_default_folder = questdlg({'Default running path not set','Select now?'},'Deafault running path','Yes', 'No','Yes');
        if strcmp(set_default_folder,'Yes')
            default_running_folder = uigetdir([],'Select Default Running Path');
            setenv('DEFAULT_RUNNING_FOLDER',default_running_folder);
            if ~isempty(default_running_folder)
                uiwait(msgbox('Successfully set default path.','Success','modal'));
            end
        end
    end
end
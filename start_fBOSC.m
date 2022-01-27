
function start_fBOSC
% Get location of the function
p = mfilename('fullpath');
folder_path = fileparts(p);

% Add relevent folders
addpath(genpath(fullfile(folder_path)));

try
    py.importlib.import_module('fooof');
    disp('Successfully imported fooof from python!');
catch
    error('Problem loading fooof from python');
end

end



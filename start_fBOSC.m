
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
    warning('You will need to use the MATLAB implementation of FOOOF');   
end

end



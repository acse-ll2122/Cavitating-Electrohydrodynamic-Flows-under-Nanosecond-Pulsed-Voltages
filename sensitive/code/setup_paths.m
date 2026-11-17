function restore=setup_paths()
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(fullfile(root,'algorithm','code'));restore=study_paths();
end

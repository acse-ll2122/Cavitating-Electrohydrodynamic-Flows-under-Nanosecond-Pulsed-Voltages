function restore=study_paths()
% Activate the single shared model; restore the caller's search path on exit.
root=project_root();
saved=path;restore=onCleanup(@()path(saved));
model=fullfile(root,'algorithm','code','model');
addpath(fullfile(root,'algorithm','code'),model,fullfile(model,'physics'),fullfile(model,'solver'),fullfile(model,'output'),'-begin');
end

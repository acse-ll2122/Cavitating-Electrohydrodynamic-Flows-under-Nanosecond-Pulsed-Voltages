function root=project_root()
% Project location independent of the caller's current folder.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
end

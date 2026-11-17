function restore = activate_code(root)
    %ACTIVATE_CODE Isolate functions from other model releases.
    savedPath=path;
    restore=onCleanup(@()path(savedPath));
    folders={'physics','solver','output'};
    names={};
    for k=1:numel(folders)
        files=dir(fullfile(root,folders{k},'*.m'));
        for j=1:numel(files)
            [~,names{end+1}]=fileparts(files(j).name);
            %#ok<AGROW>
        end
    end
    clear(names{:});
    for k=1:numel(folders),
        addpath(fullfile(root,folders{k}),'-begin');
    end
end

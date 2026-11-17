function signature = case_signature(C,root)
    %CASE_SIGNATURE Bind cached results to settings and all executable sources.
    settings=C;
    for name={'voltage','sourceDebugFile','invalidStateFile'}
        if isfield(settings,name{1}),
            settings=rmfield(settings,name{1});
        end
    end
    signature.settings=settings;
    signature.voltage_function=func2str(C.voltage);
    signature.voltage_samples=arrayfun(C.voltage,[0 .25 .5 .75 1]*1e-9);
    files=dir(fullfile(root,'**','*.m'));
    names={};
    sources={};
    for k=1:numel(files)
        file=fullfile(files(k).folder,files(k).name);
        relative=strrep(file(numel(root)+2:end),'\','/');
        if contains(relative,'/results/')||startsWith(relative,'validation/'),
            continue;
        end
        names{end+1}=relative;
        %#ok<AGROW>
        sources{end+1}=fileread(file);
        %#ok<AGROW>
    end
    [signature.names,order]=sort(names);
    signature.sources=sources(order);
end

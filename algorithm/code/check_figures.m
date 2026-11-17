function records=check_figures()
% Load without showing batch windows, retaining the saved Figure.Visible value.
root=project_root();scope=study_paths(); %#ok<NASGU>
files=dir(fullfile(root,'**','*.fig'));records=cell(numel(files),1);
for k=1:numel(files)
    file=fullfile(files(k).folder,files(k).name);[h,properties]=hgload(file,struct('Visible','off'));
    if iscell(properties),properties=properties{1};end
    visible=properties.Visible;set(h,'Units','inches');position=get(h,'Position');checked=0;
    lines=findobj(h,'Type','line');
    for j=1:numel(lines)
        data=lines(j).UserData;
        if ~isstruct(data)||~isfield(data,'source')||~isfield(data,'field'),continue;end
        name=files(k).name;
        if startsWith(data.source,'MG_epsilon/')||startsWith(data.source,'动态/')
            relative=strrep(data.source,'动态/','MG_epsilon/');relative=strrep(relative,'电压比较/','voltage/');relative=strrep(relative,'升压时间比较/','rise_time/');
            source=fullfile(root,relative);
        elseif contains(name,'epsilon')
            [~,id]=fileparts(data.source);token=regexp(id,'^(.*)_([0-9.]+)ns$','tokens','once');
            source=epsilon_profile(token{1},str2double(token{2}));
        else
            source=fullfile(root,'sensitive',strrep(data.source,'result/','results/'));
        end
        P=readtable(source);idx=1:height(P);
        if isfield(data,'indices'),idx=data.indices;end
        if ~isequal(lines(j).XData(:),P.r_um(idx))||~isequal(lines(j).YData(:),P.(data.field)(idx)),error('Check:FigureData','Saved FIG differs from current CSV: %s',file);end
        checked=checked+1;
    end
    close(h);
    if ~strcmp(visible,'on'),error('Check:FigureVisible','Saved Figure.Visible is off: %s',file);end
    records{k}=struct('file',file(numel(root)+2:end),'Visible',visible,'canvas_inches',position(3:4),'native_data_series_checked',checked);
end
out=fullfile(root,'algorithm','results','summary','verification');if ~isfolder(out),mkdir(out);end
write_json(fullfile(out,'FIG可见性核查.json'),records);
fprintf('CHECK_FIGURES_PASS %d FIG files\n',numel(files));
end

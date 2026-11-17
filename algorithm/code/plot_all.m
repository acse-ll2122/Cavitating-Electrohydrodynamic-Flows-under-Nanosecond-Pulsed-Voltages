function plot_all()
% Replot all paper figures from saved data; no PDE integration.
root=project_root();saved=pwd;cleanup=onCleanup(@()cd(saved));
scope=study_paths(); %#ok<NASGU>
T=readtable(catalog_file(),'TextType','string');
for k=1:height(T)
 S=jsondecode(fileread(fullfile(root,char(T.folder(k)),'status.json')));
 if ~S.completed,error('EHD:Incomplete','Incomplete result: %s',T.key(k));end
end
for k=1:4,plot_algorithm(k);end
cd(fullfile(root,'sensitive','code'));clear setup_paths analyze_results;
analyze_results();plot_parameters();for k=1:5,plot_groups(k);end
cd(fullfile(root,'MG_epsilon','voltage','code'));clear setup_paths analyze_results;
analyze_results();plot_epsilon('voltage');analyze_force();plot_force();
cd(fullfile(root,'MG_epsilon','rise_time','code'));clear setup_paths;
plot_epsilon('rise_time');
refresh_results();
end

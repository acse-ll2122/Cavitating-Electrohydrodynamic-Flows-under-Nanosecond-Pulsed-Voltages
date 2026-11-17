function run_all(worker,workers)
% From the project root: run_all recomputes every case, then plots all results.
% Optional worker/workers partitions the catalog into independent batches.
if nargin<1,worker=1;workers=1;end
root=project_root();scope=study_paths(); %#ok<NASGU>
T=readtable(catalog_file(),'TextType','string');
% Long voltage ramps start first so independent batches finish together.
[~,order]=sort(T.rise_ns,'descend');T=T(order,:);
for k=worker:workers:height(T),run_study_case(char(T.key(k)));end
if workers==1,check_results();plot_all();check_figures();end
end

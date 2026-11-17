function result=read_case(id,time_ns)
if nargin<1,id='BASE';end
scope=setup_paths(); %#ok<NASGU>
[C,row]=study_config(['parameter_' char(id)]);if nargin<2,time_ns=C.tFinal*1e9;end
root=project_root();folder=fullfile(root,char(row.folder));
status=jsondecode(fileread(fullfile(folder,'status.json')));
if ~status.completed,error('EHD:Incomplete','This case has not completed its current run.');end
result=load(fullfile(folder,sprintf('state_%g_ns.mat',time_ns)));result.status=status;result.source_folder=folder;
end

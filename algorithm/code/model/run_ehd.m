function status=run_ehd(model_id,U_kV,N,H_ps)
% Optional single-model run; run_all in algorithm/code reproduces the paper set.
if nargin<1,model_id='04_dynamic';end
if nargin<2,U_kV=50;end
if nargin<3,N=1600;end
if nargin<4,H_ps=.125;end
root=fileparts(mfilename('fullpath'));scope=activate_code(root); %#ok<NASGU>
[C,model]=model_config(model_id,U_kV,N,H_ps);
tag=sprintf('%s_%gkV_N%d_H%g',model_id,U_kV,N,H_ps);
out=fullfile(fileparts(fileparts(fileparts(root))),'algorithm','results','custom',model.folder,tag);
status=run_case(C,tag,out);
signature=case_signature(C,root);save(fullfile(out,'case_signature.mat'),'signature');
end

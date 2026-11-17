function status=run_ehd(model_id,U_kV,N,H_ps)
if nargin<1,model_id='04_dynamic';end
if nargin<2,U_kV=50;end
if nargin<3,N=1600;end
if nargin<4,H_ps=.125;end
scope=setup_paths(); %#ok<NASGU>
[C,~]=model_config(model_id,U_kV,N,H_ps);
root=fileparts(fileparts(mfilename('fullpath')));tag=sprintf('%s_%gkV_N%d_H%g',model_id,U_kV,N,H_ps);
status=run_case(C,tag,fullfile(root,'results','custom',tag));
end

function status=run_model(id)
if nargin<1,id='BASE';end
scope=setup_paths(); %#ok<NASGU>
if strcmp(id,'BASE')&&strcmp('epsilon','epsilon'),id='V50_T1_dynamic';end
status=run_study_case(['epsilon_' char(id)]);
end

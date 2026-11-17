function status=run_model(id)
if nargin<1,id='BASE';end
scope=setup_paths(); %#ok<NASGU>
status=run_study_case(['parameter_' char(id)]);
end

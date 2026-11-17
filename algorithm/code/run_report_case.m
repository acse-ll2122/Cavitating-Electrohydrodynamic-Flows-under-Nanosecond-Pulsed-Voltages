function status=run_report_case(id)
root=project_root(); %#ok<NASGU>
T=readtable(catalog_file(),'TextType','string');i=find(T.id==string(id)&(T.kind=="algorithm"|T.kind=="validation"));
if numel(i)~=1,error('EHD:CaseID','Choose an algorithm/validation ID from case_catalog.csv.');end
status=run_study_case(char(T.key(i)));
end

function [file,relative]=epsilon_profile(id,time_ns)
% Resolve the canonical profile through the catalog, including the shared pair.
T=readtable(catalog_file(),'TextType','string');i=find(T.kind=="epsilon"&T.id==string(id));
if numel(i)~=1,error('EHD:CaseID','Unknown dielectric comparison case.');end
relative=fullfile(fileparts(fileparts(char(T.folder(i)))),'profiles',sprintf('%s_%gns.csv',id,time_ns));
relative=strrep(relative,'\','/');file=fullfile(project_root(),relative);
end

function [C,row]=study_config(key)
% All study settings are explicit inputs; no saved result is read here.
root=project_root();T=readtable(catalog_file(),'TextType','string');
i=find(T.key==string(key));
if numel(i)~=1,error('EHD:CaseID','Choose a key from case_catalog.csv.');end
row=T(i,:);id=char(row.model_id);stage=str2double(id(1:2));
[C,~]=model_config(sprintf('%02d_fixed',stage),row.U_kV,row.N,row.H_ps);
if contains(id,'dynamic'),C.epsilonMode='maxwell-garnett';end
C.variant=['B' id];C.modelName=[C.variant ' / ' C.timeMode];
C.vofCc=row.Cc;C.vofCv=row.Cv;C.vofDNuc=row.d_nm*1e-9;C.vofN0=row.n0_m3;C.pThreshold=row.pth_MPa*1e6;
C.tRamp=row.rise_ns*1e-9;C.tFinal=C.tRamp;U=C.Umax;T=C.tRamp;C.voltage=@(t)U*(t/T);
end

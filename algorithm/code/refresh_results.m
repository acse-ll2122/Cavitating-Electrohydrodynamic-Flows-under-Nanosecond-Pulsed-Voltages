function summary=refresh_results()
% One catalog, one table of endpoint quantities and provenance.
root=project_root();T=readtable(catalog_file(),'TextType','string');
rows=cell(height(T),12);entries=cell(height(T),1);cflRows=cell(height(T),8);
for k=1:height(T)
 folder=fullfile(root,char(T.folder(k)));A=load(fullfile(folder,sprintf('state_%g_ns.mat',T.rise_ns(k))));S=jsondecode(fileread(fullfile(folder,'status.json')));
 [u,iu]=min(A.state.u);[p,ip]=min(A.state.ptotal);
 rows(k,:)={char(T.key(k)),u,A.grid.centers(iu)*1e6,p*1e-6,A.grid.centers(ip)*1e6,min(A.state.alpha),max(A.state.Er)*1e-9,S.steps,S.rejected_steps,S.min_dt_ps,S.max_relative_mass_budget_residual,A.time*1e9};
 files=dir(fullfile(folder,'state_*_ns.mat'));entries{k}=struct('key',char(T.key(k)),'id',char(T.id(k)),'kind',char(T.kind(k)),'folder',char(T.folder(k)),'endpoint',sprintf('state_%g_ns.mat',T.rise_ns(k)),'snapshots',{{files.name}},'status','completed');
 clock=load(fullfile(folder,'accepted_clock.mat'));common=~ismissing(T.clock(k))&&strlength(T.clock(k))>0;
 cflRows(k,:)={char(T.key(k)),A.C.cfl,'abs(u)+c_fr',T.H_ps(k),S.steps,S.rejected_steps,min(clock.schedule(:,2))*1e12,common};
end
summary=cell2table(rows,'VariableNames',{'key','u_min_m_s','r_u_min_um','ptotal_min_MPa','r_ptotal_min_um','alpha_min','E_max_GV_m','accepted_steps','rejected_steps','min_dt_ps','mass_budget_residual','time_ns'});
out=fullfile(root,'algorithm','results','summary');if ~isfolder(out),mkdir(out);end
writetable(summary,fullfile(out,'全部47组_当前结果指标.csv'));

writetable(T,fullfile(out,'当前47组算例索引.csv'));
cflTable=cell2table(cflRows,'VariableNames',{'key','CFL_coefficient','wave_speed','nominal_dt_ps','accepted_steps','rejected_steps','min_actual_dt_ps','common_clock_checked_each_step'});
writetable(cflTable,fullfile(out,'全部47组_CFL限步记录.csv'));
end

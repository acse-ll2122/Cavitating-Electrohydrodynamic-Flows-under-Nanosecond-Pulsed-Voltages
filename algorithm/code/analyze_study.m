function summary=analyze_study(kind)
% Endpoint statistics for one catalog category.
root=project_root();T=readtable(catalog_file(),'TextType','string');T=T(T.kind==string(kind),:);
rows=cell(height(T),12);
for k=1:height(T)
    folder=fullfile(root,char(T.folder(k)));P=readtable(fullfile(folder,sprintf('profile_%g_ns.csv',T.rise_ns(k))));
    S=jsondecode(fileread(fullfile(folder,'status.json')));[u,iu]=min(P.u);[p,ip]=min(P.ptotal_MPa);near=P.r_um<=44;
    rows(k,:)={char(T.id(k)),u,P.r_um(iu),p,P.r_um(ip),min(P.alpha),excess_tv(P.u(near)),excess_tv(P.ptotal_MPa(near)),S.steps,S.rejected_steps,S.min_dt_ps,S.max_relative_mass_budget_residual};
end
summary=cell2table(rows,'VariableNames',{'id','u_min','r_u_min_um','ptotal_min_MPa','r_ptotal_min_um','alpha_min','nearwall_u_excess_TV','nearwall_pressure_excess_TV_MPa','accepted','rejected','min_actual_dt_ps','max_relative_mass_budget_residual'});
if strcmp(kind,'parameter')
 out=fullfile(root,'sensitive','results','metrics');if ~isfolder(out),mkdir(out);end
 writetable(summary,fullfile(out,'recomputed_summary.csv'));
else
 groups={'voltage','rise_time'};masks={T.rise_ns==1,T.U_kV==50};
 for j=1:2
  out=fullfile(root,'MG_epsilon',groups{j},'results','metrics');if ~isfolder(out),mkdir(out);end
  writetable(summary(masks{j},:),fullfile(out,'recomputed_summary.csv'));
 end
end
end
function value=excess_tv(y)
value=max(0,sum(abs(diff(y)))-abs(y(end)-y(1)));
end

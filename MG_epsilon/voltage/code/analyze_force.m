function analyze_force()
% Analyze current fixed and dynamic50kV/1ns BASE force densities; no PDE update.
code=fileparts(mfilename('fullpath'));package=fileparts(code);out=fullfile(package,'results','force');if ~isfolder(out),mkdir(out);end
scope=setup_paths();root=project_root();cleanup=study_paths();records={}; %#ok<NASGU>
for mode={'fixed','dynamic'}
    name=mode{1};if strcmp(name,'fixed'),rel='sensitive/results/cases/BASE';else,rel='MG_epsilon/voltage/results/cases/V50_T1_dynamic';end;file=fullfile(root,rel,'state_1_ns.mat');B=load(file);C=B.C;g=B.grid;s=B.state;
    assert(C.Umax==50000&&C.Nr==1600&&B.time==1e-9&&C.tFinal==1e-9);
    assert(isequaln(s.pressureElectric,-.5*C.alphaEhd.*s.epsilon.*s.Er.^2));
    [f1,f2,total,check]=force_split(C,g,s);
    T=table(g.centers*1e6,f1,f2,total,f1*1e-12,f2*1e-12,total*1e-12,...
        'VariableNames',{'r_um','F1_N_m3','F2_N_m3','F_N_m3','F1_TN_m3','F2_TN_m3','F_TN_m3'});
    writetable(T,fullfile(out,[name '_F_split.csv']));
    targets=[40 41 42 45 50 60 80 100 200];rows=cell(1,numel(targets));
    for j=1:numel(targets)
        [~,i]=min(abs(g.centers*1e6-targets(j)));
        rows{j}=struct('requested_r_um',targets(j),'r_um',g.centers(i)*1e6,'F1_TN_m3',f1(i)*1e-12,'F2_TN_m3',f2(i)*1e-12,'F_TN_m3',total(i)*1e-12);
    end
    [v1,i1]=min(f1);[v2,i2]=min(f2);[vt,it]=min(total);[a1,ia]=max(abs(f1));
    records{end+1}=struct('mode',name,'source',[rel '/state_1_ns.mat'],'alphaEhd',C.alphaEhd,'epsilon0',C.epsilon0,'N',C.Nr,'time_ns',B.time*1e9,'Cc',C.vofCc,'Cv',C.vofCv,'d_nm',C.vofDNuc*1e9,'n0_m3',C.vofN0,'pth_MPa',C.pThreshold*1e-6,'rows',{rows},'F1_min_TN_m3',v1*1e-12,'F1_min_r_um',g.centers(i1)*1e6,'F2_min_TN_m3',v2*1e-12,'F2_min_r_um',g.centers(i2)*1e6,'F_min_TN_m3',vt*1e-12,'F_min_r_um',g.centers(it)*1e6,'F1_abs_max_TN_m3',a1*1e-12,'F1_abs_max_r_um',g.centers(ia)*1e6,'F1_share_cell1',f1(1)/total(1),'verification',check); %#ok<AGROW>
    if strcmp(name,'fixed'),assert(all(f1==0)&&isequaln(f2,total));end
end
fid=fopen(fullfile(out,'analysis.json'),'w');assert(fid>=0);fprintf(fid,'%s',jsonencode(struct('records',{records},'new_PDE_integrations',0,'pressure','pE=-beta*epsilon0*epsilon_r*E^2/2 in both modes','force','F1=-epsilon0*E^2*grad(epsilon_r)/2; F2=beta*epsilon0*grad(epsilon_r*E^2)/2','display_units','TN/m^3 = 10^12 N/m^3','cell_average_decomposition',true)));fclose(fid);
fprintf('Force split verified for two current BASE states; native1600-cell force sums retained.\n');
end

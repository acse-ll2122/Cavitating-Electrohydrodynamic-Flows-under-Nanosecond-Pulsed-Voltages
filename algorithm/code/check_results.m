function records=check_results()
% Check every saved flow state against the model equations and its CSV.
root=project_root();scope=study_paths(); %#ok<NASGU>
T=readtable(catalog_file(),'TextType','string');records=cell(height(T),1);
for k=1:height(T)
    [C,row]=study_config(T.key(k));folder=fullfile(root,char(row.folder));S=jsondecode(fileread(fullfile(folder,'status.json')));
    if ~S.completed||abs(S.last_time_ns-C.tFinal*1e9)>1e-10,error('Check:Incomplete','Incomplete case %s',row.key);end
    times=sscanf(char(row.snapshots_ns),'%f')';maximum=zeros(1,7);
    for t=times
        A=load(fullfile(folder,sprintf('state_%g_ns.mat',t)));s=A.state;g=A.grid;
        settings=fieldnames(C);settings(strcmp(settings,'voltage'))=[];
        for j=1:numel(settings)
            if ~isequaln(A.C.(settings{j}),C.(settings{j})),error('Check:Input','Saved input differs from the catalog: %s / %s',row.key,settings{j});end
        end
        sampleTimes=[0 .5 1]*C.tFinal;
        if ~isequal(arrayfun(A.C.voltage,sampleTimes),arrayfun(C.voltage,sampleTimes)),error('Check:Voltage','Saved voltage protocol differs: %s',row.key);end
        values=[s.alpha;s.rhoLiquid;s.rhoMixture;s.u;s.p;s.Er;s.epsilon;s.ptotal];
        if numel(g.centers)~=C.Nr||~isreal(values)||any(~isfinite(values))||any(s.alpha<=0|s.alpha>1+2*sqrt(eps))||any(s.rhoLiquid<=0)||abs(A.time-t*1e-9)>1e-21,error('Check:State','Invalid state %s at %g ns',row.key,t);end
        [~,field]=dielectric(C,g,s.alpha,A.time);rate=phase_rates(C,s.rhoLiquid,s.alpha,s.rhoMixture,s.ptotal);
        expected={s.alpha.*s.rhoLiquid+(1-s.alpha)*C.rhoVapor,liquid_pressure(s.rhoLiquid,C),-.5*C.alphaEhd.*s.epsilon.*s.Er.^2,s.p+s.pressureElectric,permittivity(s.alpha,C),field,rate};
        actual={s.rhoMixture,s.p,s.pressureElectric,s.ptotal,s.epsilon,s.Er,s.massTransferRate};
        for j=1:numel(expected),maximum(j)=max(maximum(j),max(abs(actual{j}-expected{j}))/max(realmin,max(abs(expected{j}))));end
        P=readtable(fullfile(folder,sprintf('profile_%g_ns.csv',t)));native=profile_table(C,g,s);
        difference=abs(table2array(P)-table2array(native));scale=max(1,max(abs(table2array(native)),[],1));
        if any(max(difference,[],1)./scale>1e-12),error('Check:CSV','CSV differs from MAT: %s',row.key);end
    end
    if any(maximum>1e-12)||S.max_relative_mass_budget_residual>1e-9,error('Check:Closure','Model/budget check failed: %s',row.key);end
    records{k}=struct('key',char(row.key),'snapshots',numel(times),'closure_relative_error',maximum,'mass_budget_residual',S.max_relative_mass_budget_residual,'completed',true);
end
out=fullfile(root,'algorithm','results','summary','verification');if ~isfolder(out),mkdir(out);end
write_json(fullfile(out,'全部状态与CSV核查.json'),records);
fprintf('CHECK_RESULTS_PASS %d cases\n',height(T));
end

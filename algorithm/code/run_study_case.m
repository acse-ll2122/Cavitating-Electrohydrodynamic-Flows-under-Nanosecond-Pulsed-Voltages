function status=run_study_case(key)
% Recompute one catalog entry from the pure-liquid initial condition.
root=project_root();scope=study_paths(); %#ok<NASGU>
[C,row]=study_config(key);out=fullfile(root,char(row.folder));if ~isfolder(out),mkdir(out);end
write_json(fullfile(out,'status.json'),struct('completed',false,'case_id',char(row.key),'started',datestr(now,30)));
snapshots=sscanf(char(row.snapshots_ns),'%f')'*1e-9;
if ~ismissing(row.clock)&&strlength(row.clock)>0
    schedule=table2array(readtable(fullfile(root,char(row.clock))));
    [status,trace,failures]=replay_case(C,schedule,snapshots,@snapshot,@progress);
else
    [status,trace,failures,schedule]=solver(C,char(row.key),snapshots,@snapshot,@progress);
end
save(fullfile(out,'accepted_clock.mat'),'schedule');
writetable(array2table(trace,'VariableNames',{'time_ns','dt_ps','iterations','scaled_residual','source_residual'}),fullfile(out,'time_trace.csv'));
writetable(cell2table(failures,'VariableNames',{'time_ns','dt_ps','identifier','message'}),fullfile(out,'rejections.csv'));
signature=case_signature(C,fullfile(root,'algorithm','code','model'));save(fullfile(out,'case_signature.mat'),'signature');
write_json(fullfile(out,'status.json'),status);
fprintf('FINISHED %s %.9g ns, %d steps\n',key,status.last_time_ns,status.steps);
    function snapshot(requested,grid,state,time,clockCarry,initialMixtureMass)
        save(fullfile(out,sprintf('state_%g_ns.mat',requested*1e9)),'C','grid','state','time','clockCarry','initialMixtureMass');
        P=profile_table(C,grid,state);writetable(P,fullfile(out,sprintf('profile_%g_ns.csv',requested*1e9)));
        if strcmp(row.kind,'parameter')||strcmp(row.kind,'epsilon')
            profiles=fullfile(fileparts(fileparts(out)),'profiles');if ~isfolder(profiles),mkdir(profiles);end
            writetable(P,fullfile(profiles,sprintf('%s_%gns.csv',row.id,requested*1e9)));
        end
    end
    function progress(p)
        fprintf('%s: %.6g ns, %d steps, %d rejected\n',key,p.time_ns,p.steps,p.rejected_steps);
    end
end

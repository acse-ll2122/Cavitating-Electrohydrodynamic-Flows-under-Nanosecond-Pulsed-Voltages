function status=run_case(C,tag,out)
% Save native states and explicit physical columns from a fresh integration.
if ~isfolder(out),mkdir(out);end
write_json(fullfile(out,'status.json'),struct('completed',false,'case_id',tag,'started',datestr(now,30)));
[status,trace,failures,schedule]=solver(C,tag,[0 .5 1]*C.tFinal,@snapshot,@progress);
save(fullfile(out,'accepted_clock.mat'),'schedule');
writetable(array2table(trace,'VariableNames',{'time_ns','dt_ps','iterations','scaled_residual','source_residual'}),fullfile(out,'time_trace.csv'));
writetable(cell2table(failures,'VariableNames',{'time_ns','dt_ps','identifier','message'}),fullfile(out,'rejections.csv'));
write_json(fullfile(out,'status.json'),status);
    function snapshot(requested,grid,state,time,clockCarry,initialMixtureMass)
        save(fullfile(out,sprintf('state_%g_ns.mat',requested*1e9)),'C','grid','state','time','clockCarry','initialMixtureMass');
        writetable(profile_table(C,grid,state),fullfile(out,sprintf('profile_%g_ns.csv',requested*1e9)));
    end
    function progress(p)
        fprintf('%s %.6g ns, %d steps\n',tag,p.time_ns,p.steps);
    end
end

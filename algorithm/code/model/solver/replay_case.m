function [status,trace,failures]=replay_case(C,schedule,snapshots,onSnapshot,onProgress)
% Prescribed common clocks isolate the numerical change in S1 and S07.
% These CSVs are comparison inputs, not cached flow states.
maxNumCompThreads(1);grid=setup_grid(C);
p=struct('alpha',ones(C.Nr,1),'rhoLiquid',C.rho0*ones(C.Nr,1),'u',zeros(C.Nr,1));
state=evaluate_state(C,grid,p,0);state.conservativeVaporMass=zeros(C.Nr,1);state.conservativeMomentum=zeros(C.Nr,1);
time=0;clockCarry=0;initial=sum(state.rhoMixture.*grid.volume);boundary=0;budget=0;timer=tic;k=0;
for target=snapshots
    while target-time>10*eps(max(target,realmin))
        k=k+1;h=schedule(k,2);speed=abs(state.u)+state.soundSpeed.*sqrt(state.rhoLiquid./(state.alpha.*state.rhoMixture));
        cap=C.cfl*min(grid.volume./((grid.area(1:end-1)+grid.area(2:end)).*speed));
        if h>cap*(1+32*eps),error('EHD:CFLReplay','Common clock exceeds the c_fr CFL cap.');end
        [endpoint,carry]=sum_increment(time,clockCarry,h);
        [state,d]=time_step(C,grid,state,h,endpoint);time=endpoint;clockCarry=carry;
        boundary=boundary+d.boundaryMassIncrement;budget=max(budget,abs(sum(state.rhoMixture.*grid.volume)-initial-boundary)/initial);
        if mod(k,100)==0,onProgress(struct('time_ns',time*1e9,'steps',k,'rejected_steps',0));end
    end
    onSnapshot(target,grid,state,time,clockCarry,initial);
end
status=struct('completed',true,'last_time_ns',time*1e9,'steps',k,'rejected_steps',0,'min_dt_ps',min(schedule(:,2))*1e12,'runtime_seconds',toc(timer),'max_relative_mass_budget_residual',budget,'clock','prescribed common comparison clock');
trace=zeros(0,5);failures=cell(0,4);
end

function [status,trace,failures,schedule] = solver(C,tag,snapshots,onSnapshot,onProgress)
    %SOLVER Initial state, CFL cap, time stepping and rejection diagnostics.
    % Output callbacks receive accepted states; this function writes no files.
    maxNumCompThreads(1);
    grid=setup_grid(C);
    primitive.alpha=ones(C.Nr,1);
    primitive.rhoLiquid=C.rho0*ones(C.Nr,1);
    primitive.u=zeros(C.Nr,1);
    state=evaluate_state(C,grid,primitive,0);
    state.conservativeVaporMass=zeros(C.Nr,1);
    state.conservativeMomentum=zeros(C.Nr,1);
    time=0;
    clockCarry=0;
    steps=0;
    rejected=0;
    timer=tic;
    trace=zeros(0,5);
    failures=cell(0,4);
    fieldIterations=0;
    fieldDifference=0;
    sourceMax=0;
    scaledMax=0;
    maxIterations=0;
    spaceCalls=0;
    nextReport=100;
    minDt=inf;
    initialMixtureMass=sum(state.rhoMixture.*grid.volume);
    massDriftMaxAccepted=0;
    boundaryMass=0;
    budgetResidual=0;
    schedule=zeros(0,4);
    status=struct('case_id',tag,'completed',false,'last_saved_ns',0,'error','');
        for snapshot=snapshots
            while snapshot-time>10*eps(max(snapshot,realmin))
                dt=min(C.dt,snapshot-time);
                speed=abs(state.u)+state.soundSpeed.*sqrt(state.rhoLiquid./(state.alpha.*state.rhoMixture));
                cflDt=C.cfl*min(grid.volume./((grid.area(1:end-1)+grid.area(2:end)).*speed));
                dt=min(dt,cflDt);
                [newTime,~]=sum_increment(time,clockCarry,dt);
                [candidate,ss,newCarry]=advance_step(C,grid,state,dt,newTime,clockCarry,0);
                newTime=candidate.time;
                rejected=rejected+ss.rejected;
                if size(failures,1)<256,
                    failures=[failures;ss.failures];
                    if size(failures,1)>256,
                        failures=failures(1:256,:);
                    end;
                end %#ok<AGROW>
                state=candidate;
                schedule=[schedule;ss.acceptedClock]; %#ok<AGROW>
                boundaryMass=boundaryMass+ss.boundaryMassIncrement;
                budgetResidual=max(budgetResidual,abs(sum(state.rhoMixture.*grid.volume)-initialMixtureMass-boundaryMass)/initialMixtureMass);
                massDriftMaxAccepted=max(massDriftMaxAccepted,abs((sum(state.rhoMixture.*grid.volume)-initialMixtureMass)/initialMixtureMass));
                time=newTime;
                clockCarry=newCarry;
                steps=steps+ss.accepted;
                sourceMax=max(sourceMax,ss.sourceResidual);
                scaledMax=max(scaledMax,ss.scaledResidual);
                maxIterations=max(maxIterations,ss.iterations);
                spaceCalls=spaceCalls+ss.spaceCalls;
                minDt=min(minDt,ss.minDt);
                fieldIterations=max(fieldIterations,ss.fieldIterations);
                fieldDifference=max(fieldDifference,ss.fieldMassDifference);
                if steps>=nextReport
                    nextReport=100*(floor(steps/100)+1);
                    progress=struct('time_ns',time*1e9,'steps',steps,'rejected_steps',rejected,'elapsed_s',toc(timer));
                    onProgress(progress);
                end
                if mod(steps,128)==0,
                    trace=[trace;[time*1e9,dt*1e12,ss.iterations,ss.scaledResidual,ss.sourceResidual]];
                end %#ok<AGROW>
            end
            onSnapshot(snapshot,grid,state,time,clockCarry,initialMixtureMass);
            status.last_saved_ns=snapshot*1e9;
        end
        status.completed=true;
    status.last_time_ns=time*1e9;
    status.steps=steps;
    status.rejected_steps=rejected;
    status.runtime_seconds=toc(timer);
    status.min_dt_ps=minDt*1e12;
    status.max_source_coordinate_residual=sourceMax;
    status.max_stage_scaled_residual=scaledMax;
    status.max_stage_iterations=maxIterations;
    status.space_calls=spaceCalls;
    status.max_relative_mass_drift_accepted_steps=massDriftMaxAccepted;
    status.final_relative_mass_drift=(sum(state.rhoMixture.*grid.volume)-initialMixtureMass)/initialMixtureMass;
    status.method=C.timeMode;
    status.variant=C.variant;
    status.boundary_mode=C.boundaryMode;
    status.max_field_iterations=fieldIterations;
    status.max_field_mass_difference=fieldDifference;
    status.epsilon_mode=C.epsilonMode;
    status.reconstruction=C.reconstruction;
    status.max_relative_mass_budget_residual=budgetResidual;
    status.boundary_mass_increment=boundaryMass;
end

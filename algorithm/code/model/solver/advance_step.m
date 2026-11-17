function [next,stats,newCarry]=advance_step(C,grid,state,dt,newTime,clockCarry,depth)
    %ADVANCE_STEP Bisect failed EHD trials; preserve the configured nominal step.
    try
        [next,stats]=time_step(C,grid,state,dt,newTime);
        [~,newCarry]=sum_increment(state.time,clockCarry,dt);
        stats.accepted=1;
        stats.rejected=0;
        stats.minDt=dt;
        stats.acceptedClock=[state.time dt next.time newCarry];
        stats.failures=cell(0,4);
        return;
    catch err
        if ~startsWith(err.identifier,'EHD:'),
            rethrow(err);
        end
        maxDepth=C.maxBisection;
        if depth>=maxDepth || dt/2<1e-21,
            rethrow(err);
        end
        failed={state.time*1e9,dt*1e12,err.identifier,err.message};
    end
    [mid,carryMid]=sum_increment(state.time,clockCarry,dt/2);
    [middle,s1,carryMid]=advance_step(C,grid,state,dt/2,mid,clockCarry,depth+1);
    [endpoint,~]=sum_increment(mid,carryMid,dt/2);
    [next,s2,newCarry]=advance_step(C,grid,middle,dt/2,endpoint,carryMid,depth+1);
    stats=s2;
    stats.acceptedClock=[s1.acceptedClock;s2.acceptedClock];
    stats.boundaryMassIncrement=s1.boundaryMassIncrement+s2.boundaryMassIncrement;
    stats.accepted=s1.accepted+s2.accepted;
    stats.rejected=1+s1.rejected+s2.rejected;
    stats.minDt=min(s1.minDt,s2.minDt);
    stats.iterations=max(s1.iterations,s2.iterations);
    stats.scaledResidual=max(s1.scaledResidual,s2.scaledResidual);
    stats.sourceResidual=max(s1.sourceResidual,s2.sourceResidual);
    stats.sourceIterations=s1.sourceIterations+s2.sourceIterations;
    stats.spaceCalls=s1.spaceCalls+s2.spaceCalls;
    stats.fieldIterations=max(s1.fieldIterations,s2.fieldIterations);
    stats.fieldMassDifference=max(s1.fieldMassDifference,s2.fieldMassDifference);
    stats.failures=[failed;s1.failures;s2.failures];
    if size(stats.failures,1)>256,
        stats.failures=stats.failures(1:256,:);
    end
end

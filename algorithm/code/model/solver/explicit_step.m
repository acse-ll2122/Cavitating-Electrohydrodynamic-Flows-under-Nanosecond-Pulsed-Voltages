function [next,boundaryMassIncrement] = explicit_step(C, grid, state, dt, newTime, timeMode)
    %EXPLICIT_STEP Forward Euler (B00/B01) or Step 2 SSPRK3 (B02).
    oldTime = newTime-dt;
    switch lower(timeMode)
        case 'fe'
            [next,boundaryMassIncrement] = euler_stage(C,grid,state,dt,newTime);
        case 'ssprk3'
            [stage1,b1] = euler_stage(C,grid,state,dt,newTime);
            [trial2,b2] = euler_stage(C,grid,stage1,dt,newTime+dt);
            stage2 = combine_states(C,grid,state,trial2,3/4,1/4,oldTime+dt/2);
            [trial3,b3] = euler_stage(C,grid,stage2,dt,newTime);
            boundaryMassIncrement=b1/6+b2/6+2*b3/3;
            next = combine_states(C,grid,state,trial3,1/3,2/3,newTime);
        otherwise
            error('Unknown time mode: %s',timeMode);
    end
end

function [next,boundaryMassIncrement] = euler_stage(C,grid,state,dt,time)
    op = space_operator(C,grid,state);
    boundaryMassIncrement=dt*(op.liquidMassFlux(1)+op.vaporMassFlux(1)-op.liquidMassFlux(end)-op.vaporMassFlux(end));
    rhoLiquidPartial = state.alpha.*state.rhoLiquid + dt*op.liquidMassRhs;
    rhoVaporPartial = (1-state.alpha).*C.rhoVapor + dt*op.vaporMassRhs;
    rhoMixture = state.rhoMixture + dt*op.mixtureMassRhs;
    momentumDensity = state.rhoMixture.*state.u + dt*op.momentumRhs;
    next = state_from_conservative(C,grid,rhoLiquidPartial,rhoVaporPartial, ...
        rhoMixture,momentumDensity,time);
end

function mixed = combine_states(C,grid,a,b,weightA,weightB,time)
    rhoLiquidPartial = weightA*(a.alpha.*a.rhoLiquid) ...
        + weightB*(b.alpha.*b.rhoLiquid);
    rhoVaporPartial = weightA*((1-a.alpha).*C.rhoVapor) ...
        + weightB*((1-b.alpha).*C.rhoVapor);
    rhoMixture = weightA*a.rhoMixture+weightB*b.rhoMixture;
    momentumDensity = weightA*(a.rhoMixture.*a.u)+weightB*(b.rhoMixture.*b.u);
    mixed = state_from_conservative(C,grid,rhoLiquidPartial,rhoVaporPartial, ...
        rhoMixture,momentumDensity,time);
end

function state = state_from_conservative(C,grid,rhoLiquidPartial,rhoVaporPartial, ...
        rhoMixture,momentumDensity,time)
    if any(~isfinite(rhoLiquidPartial)) || any(~isfinite(rhoVaporPartial)) ...
            || any(~isfinite(rhoMixture)) || any(~isfinite(momentumDensity))
        error('EHD:InvalidState','Non-finite conservative state.');
    end
    alpha = 1-rhoVaporPartial/C.rhoVapor;
    roundoffTolerance = 2*sqrt(eps(1));
    if any(alpha <= 0) || any(alpha > 1+roundoffTolerance)
        error('EHD:InvalidState', ...
            'Update produced alpha outside the roundoff-tolerant domain: min=%.17g, max=%.17g.', ...
        min(alpha),max(alpha));
    end
    if any(rhoLiquidPartial <= 0) || any(rhoMixture <= 0)
        error('EHD:InvalidState','Liquid partial density is non-positive.');
    end
    massClosureError = max(abs(rhoMixture-(rhoLiquidPartial+rhoVaporPartial)));
    if massClosureError > 1e-10*max(max(abs(rhoMixture)),1)
        error('EHD:MassClosure','Phase masses do not close mixture mass.');
    end
    primitive.alpha = alpha;
    primitive.rhoLiquid = rhoLiquidPartial./alpha;
    primitive.u = momentumDensity./rhoMixture;
    state = evaluate_state(C,grid,primitive,time);
end

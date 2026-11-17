function state = physics(C,primitive,time,epsilon,electric,chi)
    %PHYSICS Physical states and closures; no reconstruction or time stepping.
    %

    % Mixture continuity and momentum: Eqs. (1)--(2); Tait EOS: Eq. (3).
    % rho_m=alpha*rho_l+(1-alpha)*rho_v: Eq. (9).
    % Equal/opposite phase sources: Eqs. (10)--(11), alpha identity: Eq. (12).
    % Full Helmholtz force: Eq. (22), integrated in force_quadrature.m.
    % Local pE=-beta*epsilon*E^2/2: Eq. (23); ptotal=p+pE: Eq. (7).
    % For variable epsilon, full fH generally differs from -grad(pE).
    state = primitive;
    state.time = time;
    state.rhoMixture = state.alpha.*state.rhoLiquid+(1-state.alpha).*C.rhoVapor;
    state.p = liquid_pressure(state.rhoLiquid,C);
    dpDrhoLiquid = C.gamma*(C.p0+C.B)/C.rho0 ...
        .*(state.rhoLiquid/C.rho0).^(C.gamma-1);
    state.soundSpeed = sqrt(dpDrhoLiquid);
    state.epsilon = epsilon;
    state.Er = electric;
    state.pressureElectrostrictive = -electrostrictive_pressure(epsilon,electric,C);
    state.pressureElectric = -chi;
    state.ptotal = state.p+state.pressureElectric;
    [state.massTransferRate,condensation,vaporization,state.alphaNuc, ...
        state.bubbleRadius] = phase_rates(C,state.rhoLiquid,state.alpha, ...
        state.rhoMixture,state.ptotal);
    state.condensationSource = condensation./state.rhoLiquid;
    state.vaporizationSource = vaporization./state.rhoLiquid;
    state.alphaSource = state.massTransferRate./state.rhoLiquid;
    state.vaporizationActive = state.vaporizationSource<0;
    state.condensationActive = state.condensationSource>0;
    state.phaseChangeActive = state.vaporizationActive|state.condensationActive;
end

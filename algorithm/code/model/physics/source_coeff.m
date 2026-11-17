function [K,Kprime,rhoLiquid,pprime]=source_coeff(C,rhoMixture,rhoVaporPartial,evaporation)
    %SOURCE_COEFF Physical phase rate and analytic Jacobian at fixed full rho_m.
    % mDot=K*sign(gap)*sqrt(abs(gap)); derivative is with respect to vapor partial density.
    % rhoMixture=alpha*rhoLiquid+(1-alpha)*rhoVapor, including vapor mass.
    rhoVapor=C.rhoVapor;
    alpha=1-rhoVaporPartial/rhoVapor;
    rhoLiquid=(rhoMixture-rhoVaporPartial)./alpha;
    nuc=nucleation_state(C.vofN0,C.vofDNuc);
    z=nuc.alphaNuc+rhoVaporPartial/rhoVapor;
    radius=(3*z./(4*pi*C.vofN0*alpha)).^(1/3);
    base=3*rhoLiquid*rhoVapor./(rhoMixture.*radius).*sqrt(2./(3*rhoLiquid));
    dRhoLiquid_dVaporPartial=(rhoMixture-rhoVapor)./(rhoVapor*alpha.^2);
    baseLogPrime=.5*dRhoLiquid_dVaporPartial./rhoLiquid-(1./(rhoVapor*z)+1./(rhoVapor*alpha))/3;
    K=C.vofCc*base.*alpha.*(rhoVaporPartial/rhoVapor);
    Kprime=C.vofCc*base.*(-rhoVaporPartial/rhoVapor^2+alpha/rhoVapor+alpha.*(rhoVaporPartial/rhoVapor).*baseLogPrime);
    K(evaporation)=C.vofCv*base(evaporation).*alpha(evaporation).*z(evaporation);
    Kprime(evaporation)=C.vofCv*base(evaporation).* ...
        (-z(evaporation)/rhoVapor+alpha(evaporation)/rhoVapor ...
        +alpha(evaporation).*z(evaporation).*baseLogPrime(evaporation));
    pprime=C.gamma*(C.p0+C.B)*(rhoLiquid/C.rho0).^C.gamma./rhoLiquid.*dRhoLiquid_dVaporPartial;
end

function [rate,condensation,vaporization,alphaNuc,radius] = ...
        phase_rates(C,rhoLiquid,alpha,rhoMixture,pressureTotal)
    %PHASE_RATES Schnerr-Sauer/FMT, Eqs. (9)--(12); positive rate condenses.
    % rhoMixture=alpha*rhoLiquid+(1-alpha)*rhoVapor, without approximation.
    nuc = nucleation_state(C.vofN0,C.vofDNuc);
    alphaNuc = nuc.alphaNuc;
    radius = ((3/(4*pi*C.vofN0)).*(1+alphaNuc-alpha)./alpha).^(1/3);
    gap = pressureTotal-C.pThreshold;
    massScale = 3.*rhoLiquid.*C.rhoVapor./(rhoMixture.*radius);
    % Signed square roots give the continuous zero limit at gap=0.
    condensation = C.vofCc.*alpha.*(1-alpha).*massScale ...
        .*sqrt(2.*max(gap,0)./(3.*rhoLiquid));
    vaporization = -C.vofCv.*alpha.*(1+alphaNuc-alpha).*massScale ...
        .*sqrt(2.*max(-gap,0)./(3.*rhoLiquid));
    rate = condensation+vaporization;
end

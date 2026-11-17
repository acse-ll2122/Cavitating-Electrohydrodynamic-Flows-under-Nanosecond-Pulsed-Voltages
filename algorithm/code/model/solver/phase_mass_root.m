function [liquid,vapor,stats]=phase_mass_root(C,liquidPredictor,vaporPredictor,h,chi)
    %PHASE_MASS_ROOT Solve the original phase-mass balances in mass coordinates.
    % Inputs Lpred,Vpred [kg/m^3], h [s], chi [Pa]. No inverse EOS or gap coordinate.
    % Equal/opposite transfer eliminates one redundant unknown: L=Lpred+Vpred-V.
    total=liquidPredictor+vaporPredictor;
    if any(~isfinite([total;vaporPredictor;chi]))||any(total<=C.rhoVapor)||any(vaporPredictor>=C.rhoVapor)
        error('EHD:ImplicitInput','Invalid phase-mass source predictor.');
    end
    lo=C.qLower*ones(size(total));
    hi=C.rhoVapor*(1-eps(1))*ones(size(total));
    [left,~]=direct_rate(C,liquidPredictor,vaporPredictor,lo,chi);
    [right,~]=direct_rate(C,liquidPredictor,vaporPredictor,hi,chi);
    fl=lo-vaporPredictor+h*left;
    fh=hi-vaporPredictor+h*right;
    if any(fl>C.rootMassTolerance)||any(fh< -C.rootMassTolerance)
        error('EHD:ImplicitInput','No direct phase-mass root in retained alpha roundoff domain.');
    end
    vapor=min(hi,max(lo,vaporPredictor));
    coordinate=8*eps(C.rhoVapor);
    for iteration=1:96
        [rate,derivative]=direct_rate(C,liquidPredictor,vaporPredictor,vapor,chi);
        residual=vapor-vaporPredictor+h*rate;
        done=abs(residual)<=C.rootMassTolerance | hi-lo<=coordinate;
        if all(done)
            break;
        end
        active=~done;
        lower=active&residual<0;
        upper=active&~lower;
        lo(lower)=vapor(lower);
        hi(upper)=vapor(upper);
        candidate=vapor-residual./(1+h*derivative);
        bad=~isfinite(candidate)|candidate<=lo|candidate>=hi|candidate==vapor;
        candidate(bad)=(lo(bad)+hi(bad))/2;
        vapor(active)=candidate(active);
    end
    [rate,~]=direct_rate(C,liquidPredictor,vaporPredictor,vapor,chi);
    residual=vapor-vaporPredictor+h*rate;
    if any(abs(residual)>C.rootMassTolerance & hi-lo>coordinate)
        error('EHD:ImplicitConvergence','Direct phase-mass source did not converge.');
    end
    transfer=vaporPredictor-vapor;
    liquid=liquidPredictor+transfer;
    stats=[max(abs(residual)) iteration sum(abs(residual)>C.rootMassTolerance)];
end

function [rate,derivative]=direct_rate(C,liquidPredictor,vaporPredictor,vapor,chi)
    liquid=liquidPredictor+(vaporPredictor-vapor);
    alpha=1-vapor/C.rhoVapor;
    rhoLiquid=liquid./alpha;
    total=liquid+vapor;
    pressureTotal=liquid_pressure(rhoLiquid,C)-chi;
    rate=phase_rates(C,rhoLiquid,alpha,total,pressureTotal);
    % The existing analytic mass-coordinate Jacobian changes no physical law.
    gap=pressureTotal-C.pThreshold;
    [K,Kprime,~,pprime]=source_coeff(C,total,vapor,gap<0);
    derivative=Kprime.*sign(gap).*sqrt(abs(gap));
    active=gap~=0;
    derivative(active)=derivative(active)+K(active).*pprime(active)./(2*sqrt(abs(gap(active))));
    derivative(~active & K~=0)=Inf;
end

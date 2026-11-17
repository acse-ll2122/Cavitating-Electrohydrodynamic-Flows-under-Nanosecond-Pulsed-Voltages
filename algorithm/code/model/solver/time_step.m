function [next,stats]=time_step(C,g,s,h,newTime)
    %TIME_STEP Step 2 explicit stages or Step 3 IMEX-SSP3(4,3,3).
    stats=struct('iterations',0,'scaledResidual',0,'sourceResidual',0,'sourceIterations',0,'spaceCalls',1,'fieldIterations',0,'fieldMassDifference',0);
    stats.boundaryMassIncrement=0;
    if C.stage<=2
        [next,stats.boundaryMassIncrement]=explicit_step(C,g,s,h,newTime,C.timeMode);
        next.conservativeLiquidMass=next.alpha.*next.rhoLiquid;
        next.conservativeVaporMass=(1-next.alpha)*C.rhoVapor;
        next.conservativeMomentum=next.rhoMixture.*next.u;
        return;
    end
    N=C.Nr;
    rhoLiquidPartial0=s.alpha.*s.rhoLiquid;
    if isfield(s,'conservativeLiquidMass')
        rhoLiquidPartial0=s.conservativeLiquidMass;
    end
    momentumDensity0=s.rhoMixture.*s.u;
    rhoVaporPartial0=(1-s.alpha)*C.rhoVapor;
    if isfield(s,'conservativeVaporMass')
        rhoVaporPartial0=s.conservativeVaporMass;
    end
    liquidTransportRhs=zeros(N,4);
    momentumRhs=liquidTransportRhs;
    vaporTransportRhs=liquidTransportRhs;
    source=liquidTransportRhs;
    boundaryRates=zeros(1,4);
    stats.spaceCalls=3;

    % ====== Four implicit phase-source stages; explicit transport and force ======
    for i=1:4
        momentumDensity=momentumDensity0+h*momentumRhs*C.Ae(i,:)';
        rhoLiquidPredictor=rhoLiquidPartial0+h*liquidTransportRhs*C.Ae(i,:)'+h*source*C.Ai(i,:)';
        rhoVaporPredictor=rhoVaporPartial0+h*vaporTransportRhs*C.Ae(i,:)'-h*source*C.Ai(i,:)';
        ts=s.time+C.ci(i)*h;
        hi=h*C.Ai(i,i);
        [rhoLiquidPartial,rhoVaporPartial,d]=coupled_source(C,g,rhoLiquidPredictor,rhoVaporPredictor,hi,ts);
        stats.fieldIterations=max(stats.fieldIterations,d(3));
        stats.fieldMassDifference=max(stats.fieldMassDifference,d(4));
        stats.sourceResidual=max(stats.sourceResidual,d(1));
        stats.sourceIterations=stats.sourceIterations+d(2);
        source(:,i)=(rhoVaporPredictor-rhoVaporPartial)/hi;
        stage=from_mass(C,g,rhoLiquidPartial,momentumDensity,rhoVaporPartial,ts);
        if i>=2
            op=space_operator(C,g,stage);
            boundaryRates(i)=op.liquidMassFlux(1)+op.vaporMassFlux(1)-op.liquidMassFlux(end)-op.vaporMassFlux(end);
            liquidTransportRhs(:,i)=-(op.liquidMassFlux(2:end)-op.liquidMassFlux(1:end-1))./g.volume;
            momentumRhs(:,i)=op.momentumRhs;
            vaporTransportRhs(:,i)=-(op.vaporMassFlux(2:end)-op.vaporMassFlux(1:end-1))./g.volume;
        end
    end
    rhoLiquidPartial=rhoLiquidPartial0+h*(liquidTransportRhs+source)*C.weights;
    momentumDensity=momentumDensity0+h*momentumRhs*C.weights;
    rhoVaporPartial=rhoVaporPartial0+h*(vaporTransportRhs-source)*C.weights;
    stats.boundaryMassIncrement=h*boundaryRates*C.weights;
    next=from_mass(C,g,rhoLiquidPartial,momentumDensity,rhoVaporPartial,newTime);
end

function s=from_mass(C,g,rhoLiquidPartial,momentumDensity,rhoVaporPartial,t)
    rhoMixture=rhoLiquidPartial+rhoVaporPartial;
    if any(~isfinite([rhoMixture;momentumDensity;rhoVaporPartial]))||any(rhoMixture<=C.rhoVapor)||any(rhoVaporPartial<C.qLower)||any(rhoVaporPartial>=C.rhoVapor)
        error('EHD:InvalidState','IMEX stage outside alpha roundoff domain.');
    end
    p.alpha=1-rhoVaporPartial/C.rhoVapor;
    p.rhoLiquid=rhoLiquidPartial./p.alpha;
    p.u=momentumDensity./rhoMixture;
    s=evaluate_state(C,g,p,t);
    s.conservativeLiquidMass=rhoLiquidPartial;
    s.conservativeVaporMass=rhoVaporPartial;
    s.conservativeMomentum=s.rhoMixture.*s.u;
end

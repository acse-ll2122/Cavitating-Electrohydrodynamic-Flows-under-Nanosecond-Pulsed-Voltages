function [liquid,rhoVaporPartial,stats]=coupled_source(C,g,liquidPredictor,rhoVaporPredictor,h,t)
    %COUPLED_SOURCE L=Lpred+h*mDot, V=Vpred-h*mDot; equal/opposite transfer.
    guess=min(C.rhoVapor*(1-eps(1)),max(C.qLower,rhoVaporPredictor));
    maximumResidual=0;
    innerIterations=0;
    for iteration=1:C.fieldMaxIterations
        alpha=1-guess/C.rhoVapor;
        [~,~,chi]=dielectric(C,g,alpha,t);
        [liquid,rhoVaporPartial,d]=phase_mass_root(C,liquidPredictor,rhoVaporPredictor,h,chi);
        maximumResidual=max(maximumResidual,d(1));
        innerIterations=innerIterations+d(2);
        difference=max(abs(rhoVaporPartial-guess));
        if difference<=C.fieldMassTolerance
            stats=[maximumResidual innerIterations iteration difference];
            return;
        end
        guess=rhoVaporPartial;
    end
    error('EHD:FieldConvergence','Dielectric/source coupling failed: delta q = %.17g.',difference);
end

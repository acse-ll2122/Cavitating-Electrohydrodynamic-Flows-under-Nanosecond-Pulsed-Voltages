function state = evaluate_state(C,grid,primitive,time)
    %EVALUATE_STATE Retained stage-domain checks, then physical closures.
    roundoffTolerance=2*sqrt(eps(1));
    if any(primitive.alpha<=0)||any(primitive.alpha>1+roundoffTolerance)
        error('EHD:InvalidState','Alpha is outside the roundoff-tolerant model domain.');
    end
    if any(primitive.rhoLiquid<=0)
        error('EHD:InvalidState','Liquid density is non-positive.');
    end
    [epsilon,electric,chi]=dielectric(C,grid,primitive.alpha,time);
    state=physics(C,primitive,time,epsilon,electric,chi);
end

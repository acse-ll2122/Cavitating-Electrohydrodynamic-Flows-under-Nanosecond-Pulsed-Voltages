function epsilon = permittivity(alpha,C)
    %PERMITTIVITY Fixed epsilon or Maxwell-Garnett, Eq. (15); epsilon in F/m.
    if strcmpi(C.epsilonMode,'constant')
        epsilon = C.epsilonLiquid+zeros(size(alpha));
    elseif strcmpi(C.epsilonMode,'maxwell-garnett')
        f = 1-alpha;
        el = C.epsilonLiquid;
        ev = C.epsilonVapor;
        epsilon = el+3*f*el*(ev-el)./(ev+2*el-f*(ev-el));
    else
        error('EHD:EpsilonMode','Use constant or maxwell-garnett.');
    end
end

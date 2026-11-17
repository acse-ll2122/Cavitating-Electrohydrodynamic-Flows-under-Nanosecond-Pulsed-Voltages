function pressure = liquid_pressure(rhoLiquid,C)
    %LIQUID_PRESSURE Tait EOS, Model Description Eq. (4).
    % p=(p0+B)*(rho_l/rho0)^gamma-B [Pa].
    pressure = (C.p0+C.B).*(rhoLiquid./C.rho0).^C.gamma-C.B;
end

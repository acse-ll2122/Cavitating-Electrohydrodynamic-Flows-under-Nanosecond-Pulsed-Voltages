function [epsilon,electric,chi] = dielectric(C,grid,alpha,time)
    %DIELECTRIC Cellwise-epsilon voltage integral, Eq. (21).
    epsilon = permittivity(alpha,C);
    voltageIntegral = sum((1./epsilon).*(1./grid.left-1./grid.right));
    electric = electric_field(C.voltage(time),epsilon,grid.centers,voltageIntegral);
    if nargout>2
        chi = electrostrictive_pressure(epsilon,electric,C);
    end
end

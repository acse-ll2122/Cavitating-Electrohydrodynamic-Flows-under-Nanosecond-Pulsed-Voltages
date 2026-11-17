function electric = electric_field(voltage,epsilon,radius,voltageIntegral)
    %ELECTRIC_FIELD Gauss/voltage solution, Eq. (16); electric in V/m.
    % voltageIntegral = integral_a^b dr/(epsilon*r^2); quadrature is in dielectric.
    electric = (voltage/voltageIntegral)./(epsilon.*radius.^2);
end

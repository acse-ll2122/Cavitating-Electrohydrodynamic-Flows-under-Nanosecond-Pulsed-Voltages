function chi = electrostrictive_pressure(epsilon,electric,C)
    %ELECTROSTRICTIVE_PRESSURE chi=alphaEhd*epsilon*E^2/2 [Pa], Eq. (19).
    chi = .5*C.alphaEhd.*epsilon.*electric.^2;
end

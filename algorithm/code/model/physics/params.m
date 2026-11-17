function C = params()
    %PARAMS Physical constants in SI units; numerical choices are in model_config.

    % ====== Geometry and reference liquid ======
    C.a = 40e-6;
    C.b = 200e-6;
    C.waterTemperature = 293.15;
    C.rho0 = 1000;
    C.rhoVapor = .02;
    C.p0 = 1e5;
    C.B = 3.07e8;
    C.gamma = 7.5;

    % ====== Dielectric and electrostrictive closure ======
    C.epsilon0 = 8.854e-12;
    C.epsilonLiquid = 81*C.epsilon0;
    C.epsilonVapor = C.epsilon0;
    C.alphaEhd = 1.5;

    % ====== Finite mass transfer ======
    C.pThreshold = -30e6;
    C.vofN0 = 1e+20;
    C.vofDNuc = 12*1e-9;
    C.vofCc = 2;
    C.vofCv = 0.1;
end

function status=run_model(U_kV,N,H_ps)
    %RUN_MODEL 00_fixed: 0--1 ns linear ramp; shared equations and solver.
    if nargin<1,
        U_kV=50;
    end
    if nargin<2,
        N=1600;
    end
    if nargin<3,
        H_ps=.125;
    end
    root=fileparts(fileparts(mfilename('fullpath')));
    savedPath=path;
    restore=onCleanup(@()path(savedPath));
    addpath(root,'-begin');
    status=run_ehd('00_fixed',U_kV,N,H_ps);
end

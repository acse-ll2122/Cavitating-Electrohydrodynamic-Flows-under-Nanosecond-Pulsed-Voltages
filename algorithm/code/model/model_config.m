function [C,model]=model_config(model_id,U_kV,N,H_ps)
    %MODEL_CONFIG Selects algorithms; physical constants and laws stay in physics/.
    if nargin<1,
        model_id='04_dynamic';
    end
    if nargin<2,
        U_kV=50;
    end
    if nargin<3,
        N=1600;
    end
    if nargin<4,
        H_ps=.125;
    end
    if isnumeric(model_id),
        model_id=sprintf('%02d_fixed',model_id);
    end
    model_id=lower(char(model_id));
    if any(strcmp(model_id,{'00','01','02','03','04'})),
        model_id=[model_id '_fixed'];
    end
    ids={'00_fixed','01_fixed','02_fixed','03_fixed','04_fixed','04_dynamic'};
    folders={'B00_fixed','B01_fixed','B02_fixed','B03_fixed','B04_fixed','B04_dynamic'};
    ix=find(strcmp(model_id,ids),1);
    assert(~isempty(ix),'Use 00_fixed..04_fixed or 04_dynamic.');
    assert(isscalar(U_kV)&&isfinite(U_kV)&&U_kV>0);
    assert(isscalar(N)&&N==round(N)&&N>=4);
    assert(isscalar(H_ps)&&isfinite(H_ps)&&H_ps>0);
    C=params();
    C.meshBeta=.2;
    C.cfl=.35;
    C.Nr=N;
    C.dt=H_ps*1e-12;
    C.tRamp=1e-9;
    C.tFinal=1e-9;
    C.Umax=U_kV*1000;
    U_peak=C.Umax;
    ramp_time=C.tRamp;
    C.voltage=@(t)U_peak*(t/ramp_time);
    stages=[0 0 2 3 3 3];
    C.stage=stages(ix);
    C.spaceMode='quadratic';
    C.boundaryMode='quadratic';
    C.timeMode='fe';
    if ix==1,
        C.spaceMode='central';
        C.boundaryMode='cell';
    end
    if ix==3,
        C.timeMode='ssprk3';
    end
    if ix>=4,
        C.timeMode='imex_ssp3_433';
    end
    if ix>=5,
        C.boundaryMode='walls4';
    end
    C.variant=folders{ix};
    C.modelName=[folders{ix} ' / ' C.timeMode];
    C.rootMassTolerance=64*eps(C.rhoVapor);
    C.qLower=-C.rhoVapor*2*sqrt(eps(1));
    C.maxBisection=12;
    if ix==1,
        C.maxBisection=30;
    end
    C.epsilonMode='constant';
    if ix==6,
        C.epsilonMode='maxwell-garnett';
    end
    C.fieldMaxIterations=40;
    C.fieldMassTolerance=4*C.rootMassTolerance;
    a=.24169426078821;
    b=.06042356519705;
    e=.12915286960590;
    C.Ae=[0 0 0 0;0 0 0 0;0 1 0 0;0 1/4 1/4 0];
    C.Ai=[a 0 0 0;-a a 0 0;0 1-a a 0;b e .5-b-e-a a];
    C.weights=[0;1/6;1/6;2/3];
    C.ci=sum(C.Ai,2);
    C.reconstruction='quadratic';
    if ix==1,
        C.reconstruction='linear';
    end
    model=struct('id',ids{ix},'folder',folders{ix},'epsilon_mode',C.epsilonMode,'stage',C.stage);
end

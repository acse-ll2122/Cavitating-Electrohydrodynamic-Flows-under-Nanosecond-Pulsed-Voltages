function op = quadratic_operator(C,grid,state)
    %QUADRATIC_OPERATOR Steps 1 and 4; fixed three-cell reconstruction, no WENO.
    % Face bounds are retained; they are not post-step clipping of alpha.

    % ====== Conservative state: both phase densities retained ======
    N = C.Nr;
    area = grid.area;
    rhoLiquidPartial = state.alpha.*state.rhoLiquid;
    rhoVaporPartial = (1-state.alpha).*C.rhoVapor;
    rhoMixture = rhoLiquidPartial+rhoVaporPartial;
    momentumDensity = rhoMixture.*state.u;
    [WL, WR, PL, PM, PR] = quadratic_traces(grid);

    % ====== Step 1: quadratic face traces ======
    rhoLiquidPartialLeft = WL*rhoLiquidPartial;
    rhoLiquidPartialRight = WR*rhoLiquidPartial;
    rhoVaporPartialLeft = WL*rhoVaporPartial;
    rhoVaporPartialRight = WR*rhoVaporPartial;
    momentumLeft = WL*momentumDensity;
    momentumRight = WR*momentumDensity;
    aL = WL*state.alpha;
    aR = WR*state.alpha;
    uL = momentumLeft./(rhoLiquidPartialLeft+rhoVaporPartialLeft);
    uR = momentumRight./(rhoLiquidPartialRight+rhoVaporPartialRight);
    % Only velocity conditions change; density/alpha/pressure traces stay original.
    uL(1)=0;uR(1)=0;
    uL(end)=outer_velocity(state.u,grid);uR(end)=uL(end);
    momentumLeft([1 end])=(rhoLiquidPartialLeft([1 end])+rhoVaporPartialLeft([1 end])).*uL([1 end]);
    momentumRight([1 end])=(rhoLiquidPartialRight([1 end])+rhoVaporPartialRight([1 end])).*uR([1 end]);
    pL = liquid_pressure(rhoLiquidPartialLeft./aL,C);
    pR = liquid_pressure(rhoLiquidPartialRight./aR,C);
    electricAlphaFace=0.5*(aL+aR);

    % ====== Step 4: four-cell wall pressure ======
    if strcmp(C.boundaryMode,'walls4')
        pL(1) = liquid_pressure((grid.innerWeights4*rhoLiquidPartial(1:4))/(grid.innerWeights4*state.alpha(1:4)),C);
        pR(end) = liquid_pressure((grid.outerWeights4*rhoLiquidPartial(end-3:end))/(grid.outerWeights4*state.alpha(end-3:end)),C);
    end

    % ====== Rusanov mixture mass and momentum fluxes ======
    speed = abs(state.u)+state.soundSpeed.*sqrt(state.rhoLiquid./(state.alpha.*rhoMixture));
    s = zeros(N+1,1);
    s(2:N) = max(speed(1:N-1),speed(2:N));
    rhoMixtureLeft = rhoLiquidPartialLeft+rhoVaporPartialLeft;
    rhoMixtureRight = rhoLiquidPartialRight+rhoVaporPartialRight;
    mixtureDensityFlux = 0.5*(uL.*rhoMixtureLeft+uR.*rhoMixtureRight) - 0.5*s.*(rhoMixtureRight-rhoMixtureLeft);
    momentumDensityFlux = 0.5*(momentumLeft.*uL+pL+momentumRight.*uR+pR) - 0.5*s.*(momentumRight-momentumLeft);

    % ====== Vapor mass fraction: adjacent-cell and donor bounds ======
    Y = rhoVaporPartial./rhoMixture;
    YL = rhoVaporPartialLeft./(rhoLiquidPartialLeft+rhoVaporPartialLeft);
    YR = rhoVaporPartialRight./(rhoLiquidPartialRight+rhoVaporPartialRight);
    lo = [Y(1);min(Y(1:N-1),Y(2:N));Y(N)];
    hi = [Y(1);max(Y(1:N-1),Y(2:N));Y(N)];
    Yface = YL;
    Yface(mixtureDensityFlux<0) = YR(mixtureDensityFlux<0);
    Yface = min(max(Yface,lo),hi);
    Ydonor = [Y(1);Y];
    Yright = [Y;Y(N)];
    Ydonor(mixtureDensityFlux<0) = Yright(mixtureDensityFlux<0);
    Yface = min(max(Yface,0),2*max(Ydonor,0));
    vaporDensityFlux = mixtureDensityFlux.*Yface;
    liquidDensityFlux = mixtureDensityFlux-vaporDensityFlux;

    % ====== Inner no-penetration and outer Neumann velocity ======
    liquidDensityFlux(1) = 0;
    vaporDensityFlux(1) = 0;
    momentumDensityFlux(1) = pL(1);
    % Retain the complete original p3/p4 pressure; add open-boundary advection.
    momentumDensityFlux(end) = (rhoLiquidPartialRight(end)+rhoVaporPartialRight(end))*uR(end)^2+pR(end);
    liquidMassFlux = area.*liquidDensityFlux;
    vaporMassFlux = area.*vaporDensityFlux;
    momentumFlux = area.*momentumDensityFlux;
    alphaFlux = area.*0.5.*(uL+uR).*0.5.*(aL+aR);
    alphaFlux(1) = 0;

    % ====== Simpson pressure geometry and complete electric force ======
    left = grid.left;
    right = grid.right;
    mid = 0.5*(left+right);
    pAt = @(P) liquid_pressure((P*rhoLiquidPartial)./(P*state.alpha),C);
    pressureGeometry = (grid.width/6).*(2*pAt(PL).*left ...
        + 8*pAt(PM).*mid + 2*pAt(PR).*right)./grid.volume;
    [electricForce,electricalDetails] = force_quadrature(C,grid,state,electricAlphaFace,PL,PM,PR);
    massTransfer = state.massTransferRate;
    op.liquidMassRhs = -flux_divergence(grid,liquidMassFlux)+massTransfer;
    op.vaporMassRhs = -flux_divergence(grid,vaporMassFlux)-massTransfer;
    op.mixtureMassRhs = op.liquidMassRhs+op.vaporMassRhs;
    op.alphaRhs = -op.vaporMassRhs/C.rhoVapor;
    op.momentumRhs = -flux_divergence(grid,momentumFlux)+pressureGeometry+electricForce;
    op.momentumFlux = momentumFlux;
    op.pressureGeometry = pressureGeometry;
    op.electricForce = electricForce;
    op.massTransferRate = massTransfer;
    if ~strcmpi(C.epsilonMode,'constant')
        op.electricStressFace=electricalDetails.electricStressFace;
        op.electricCellIntegral=electricalDetails.electricCellIntegral;
        op.displacementSquared=electricalDetails.displacementSquared;
    end
    op.alphaFlux = alphaFlux;
    op.liquidMassFlux = liquidMassFlux;
    op.vaporMassFlux = vaporMassFlux;
end

function divF = flux_divergence(grid,flux)
    divF = (flux(2:end)-flux(1:end-1))./grid.volume;
end

function ub=outer_velocity(u,g)
% Neumann velocity fit only; no change to any density/alpha/pressure polynomial.
persistent key weights
N=numel(u);newKey=[N,g.faces(end),sum(g.faces)];
if isempty(key)||any(abs(key-newKey)>1e-20)
    st=N-2:N;face=g.faces(end);h=max(g.width(st));powers=[0 2 3];M=zeros(3);
    for k=1:3
        lo=(g.left(st(k))-face)/h;hi=(g.right(st(k))-face)/h;
        for j=1:3
            p=powers(j);M(k,j)=face^2*(hi^(p+1)-lo^(p+1))/(p+1)+2*face*h*(hi^(p+2)-lo^(p+2))/(p+2)+h^2*(hi^(p+3)-lo^(p+3))/(p+3);
        end
        M(k,:)=M(k,:)/M(k,1);
    end
    weights=(M'\[1;0;0])';key=newKey;
end
ub=weights*u(end-2:end);
end

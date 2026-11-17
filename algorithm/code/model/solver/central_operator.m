function op = central_operator(C,grid,state)
    %CENTRAL_OPERATOR B00 baseline: distance-weighted faces and central electric gradient.
    % MG baseline includes the complete product-gradient force.
    area = grid.area;
    rhoLiquidPartial = state.alpha.*state.rhoLiquid;
    rhoVaporPartial = (1-state.alpha).*C.rhoVapor;
    rhoMixture = rhoLiquidPartial+rhoVaporPartial;
    momentumDensity = rhoMixture.*state.u;
    [aL, aR] = face_states_central(grid,state.alpha);
    [lL, lR] = face_states_central(grid,rhoLiquidPartial);
    [vL, vR] = face_states_central(grid,rhoVaporPartial);
    [momL, momR] = face_states_central(grid,momentumDensity);
    [pL, pR] = face_states_central(grid,state.p);
    [gL, gR] = face_states_central(grid,grid.centers.^2.*state.u);
    gradE2 = gradient_o2(grid.centers, state.Er.^2);
    velocityFlux = 0.5*(gL+gR);
    uFace = velocityFlux./area;
    uFace(1) = 0;
    uFace(end) = state.u(end); % Constant velocity extension: du/dr=0 at b.
    alphaFace = 0.5*(aL+aR);
    alphaFlux = area .* uFace .* alphaFace;
    alphaFlux(1) = 0;
    pressureFace = 0.5*(pL+pR);
    liquidMassDensityFlux = uFace.*(0.5*(lL+lR));
    vaporMassDensityFlux = uFace.*(0.5*(vL+vR));
    momentumDensityFlux = uFace.*(0.5*(momL+momR))+pressureFace;
    liquidMassFlux = area.*liquidMassDensityFlux;
    vaporMassFlux = area.*vaporMassDensityFlux;
    liquidMassFlux(1) = 0;
    vaporMassFlux(1) = 0;
    momentumFlux = area.*momentumDensityFlux;
    momentumFlux(1) = area(1).*pressureFace(1);
    massTransfer = state.massTransferRate;
    pressureGeometry = state.p.*flux_divergence(grid,area);
    electricForce = 0.5*C.alphaEhd.*state.epsilon.*gradE2;
    if ~strcmpi(C.epsilonMode,'constant')
        % Same complete Helmholtz law as B01--B04, centrally differentiated.
        electricForce = -0.5*state.Er.^2.*gradient_o2(grid.centers,state.epsilon) ...
            +0.5*C.alphaEhd*gradient_o2(grid.centers,state.epsilon.*state.Er.^2);
    end
    op.liquidMassRhs = -flux_divergence(grid,liquidMassFlux)+massTransfer;
    op.vaporMassRhs = -flux_divergence(grid,vaporMassFlux)-massTransfer;
    op.mixtureMassRhs = op.liquidMassRhs+op.vaporMassRhs;
    op.alphaRhs = -op.vaporMassRhs/C.rhoVapor;
    op.momentumRhs = -flux_divergence(grid,momentumFlux) ...
        + pressureGeometry+electricForce;
    op.alphaFlux = alphaFlux;
    op.liquidMassFlux = liquidMassFlux;
    op.vaporMassFlux = vaporMassFlux;
end

function divF = flux_divergence(grid,flux)
    divF = (flux(2:end)-flux(1:end-1))./grid.volume;
end

function [left,right] = face_states_central(grid,cellValues)
    distance = grid.centers(2:end)-grid.centers(1:end-1);
    rightWeight = (grid.faces(2:end-1)-grid.centers(1:end-1))./distance;
    face = [cellValues(1);(1-rightWeight).*cellValues(1:end-1)+rightWeight.*cellValues(2:end);cellValues(end)];
    left = face;
    right = face;
end

function grad = gradient_o2(r,cellValues)
    grad = zeros(size(cellValues));
    grad(2:end-1) = (cellValues(3:end)-cellValues(1:end-2))./(r(3:end)-r(1:end-2));
    grad(1) = (cellValues(2)-cellValues(1))/(r(2)-r(1));
    grad(end) = (cellValues(end)-cellValues(end-1))/(r(end)-r(end-1));
end

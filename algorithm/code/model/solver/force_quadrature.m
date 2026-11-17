function [electricForce,details] = force_quadrature(C,grid,state,electricAlphaFace,PL,PM,PR)
    %FORCE_QUADRATURE Cell average of the complete Helmholtz law, Eq. (32).
    % fH=-(epsilon0/2)*E^2*grad(epsilon_r)
    %    +(alphaEhd*epsilon0/2)*grad(epsilon_r*E^2).
    % With D=epsilon*E*r^2, Gauss law gives
    % <fH>={[(alphaEhd+1)*D^2/(2*epsilon*r^2)]_faces
    %       -(alphaEhd-1)*D^2*integral dr/(epsilon*r^3)}/volume.
    % This identity integrates both terms; it is not a different physical law.
    left=grid.left;
    right=grid.right;
    mid=.5*(left+right);
    details=struct();
    if strcmpi(C.epsilonMode,'constant')
        fieldMoment = (state.Er.*grid.centers.^2).^2;
        electricForce = 0.5*C.alphaEhd.*state.epsilon ...
            .*2.*fieldMoment.*(right.^-2-left.^-2)./grid.volume;
    else
        a = C.alphaEhd;
        D2 = (state.Er(1)*state.epsilon(1)*grid.centers(1)^2)^2;
        alphaFace = min(max(electricAlphaFace,0),1);
        aAt = @(P) min(max(P*state.alpha,0),1);
        faceTerm = 0.5*(a+1)*D2./(permittivity(alphaFace,C).*grid.faces.^2);
        integral = (grid.width/6).*(1./(permittivity(aAt(PL),C).*left.^3) ...
            + 4./(permittivity(aAt(PM),C).*mid.^3) + 1./(permittivity(aAt(PR),C).*right.^3));
        electricForce = (faceTerm(2:end)-faceTerm(1:end-1)-(a-1)*D2*integral)./grid.volume;
    end
    if ~strcmpi(C.epsilonMode,'constant')
        details.electricStressFace=faceTerm;
        details.electricCellIntegral=integral;
        details.displacementSquared=D2;
    end
end

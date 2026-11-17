function [first,second,total,verification]=force_split(C,g,s)
% FORCE_SPLIT Native spherical-volume averages of the two Helmholtz terms.
% Inputs: configuration, grid and accepted state. Output units: N/m^3.
[WL,WR,PL,PM,PR]=quadratic_traces(g);alphaFace=.5*(WL*s.alpha+WR*s.alpha);
[total,detail]=force_quadrature(C,g,s,alphaFace,PL,PM,PR);
if strcmpi(C.epsilonMode,'constant')
    first=zeros(size(total));second=total;
else
    % A=[1/(epsilon_abs*r^2)]_faces, J=integral dr/(epsilon_abs*r^3).
    % Integral identities: F1=D^2*(A+2J)/(2V), F2=beta*D^2*(A-2J)/(2V).
    faceEpsilon=permittivity(min(max(alphaFace,0),1),C);
    A=diff(1./(faceEpsilon.*g.faces.^2));J=detail.electricCellIntegral;D2=detail.displacementSquared;
    first=.5*D2*(A+2*J)./g.volume;
    second=.5*C.alphaEhd*D2*(A-2*J)./g.volume;
end
error=first+second-total;scale=max(1,max(abs(total)));
assert(all(isfinite([first;second;total]))&&max(abs(error))/scale<5e-13);
verification=struct('sum_max_error_N_m3',max(abs(error)),'sum_scaled_error',max(abs(error))/scale,'units','N/m^3','cell_average',true,'reconstruction','same face epsilon and Simpson quadrature as force_quadrature');
end

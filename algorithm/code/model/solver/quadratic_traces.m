function [WL,WR,PL,PM,PR] = quadratic_traces(grid)
    %QUADRATIC_TRACES Three-cell spherical moment fits, Eqs. (24)--(25).
    persistent key cache
    N = numel(grid.centers);
    newKey = [N,grid.faces(1),grid.faces(end),sum(grid.faces)];
    if isempty(key) || any(abs(key-newKey) > 1e-20)
        stencil = @(i) max(1,min(N-2,i-1))+(0:2);
        WL = zeros(N+1,N);
        WR = zeros(N+1,N);
        PL = zeros(N,N);
        PM = zeros(N,N);
        PR = zeros(N,N);
        for i = 1:N
            st = stencil(i);
            PL(i,st) = average_weights(grid,st,grid.left(i));
            PM(i,st) = average_weights(grid,st,0.5*(grid.left(i)+grid.right(i)));
            PR(i,st) = average_weights(grid,st,grid.right(i));
        end
        WL(2:N+1,:) = PR;
        WR(1:N,:) = PL;
        WL(1,:) = PL(1,:);
        WR(N+1,:) = PR(N,:);
        cache = {sparse(WL),sparse(WR),sparse(PL),sparse(PM),sparse(PR)};
        key = newKey;
    end
    [WL,WR,PL,PM,PR] = cache{:};
end

function weights = average_weights(grid,stencil,x)
    h = max(grid.width(stencil));
    M = zeros(3,3);
    for k = 1:3
        a = (grid.left(stencil(k))-x)/h;
        b = (grid.right(stencil(k))-x)/h;
        for power = 0:2
            M(k,power+1) = x^2*(b^(power+1)-a^(power+1))/(power+1) ...
                +2*x*h*(b^(power+2)-a^(power+2))/(power+2) ...
                +h^2*(b^(power+3)-a^(power+3))/(power+3);
        end
        M(k,:) = M(k,:)/M(k,1);
    end
    weights = (M'\[1;0;0])';
end

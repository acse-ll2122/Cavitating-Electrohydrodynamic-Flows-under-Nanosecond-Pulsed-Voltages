function grid=setup_grid(C)
    %SETUP_GRID Cache geometry-only quantities without changing the power mesh.
    grid=mesh1d(C);
    % Boundary-only moment fits; internal reconstruction is unchanged.
    grid.innerWeights4=boundary_weights(grid,1:4,grid.faces(1));
    grid.outerWeights4=boundary_weights(grid,C.Nr-3:C.Nr,grid.faces(end));
end

function weights=boundary_weights(grid,stencil,face)
    % Exact spherical volume moments in a dimensionless local coordinate.
    n=numel(stencil);
    h=max(grid.width(stencil));
    M=zeros(n,n);
    for k=1:n
        a=(grid.left(stencil(k))-face)/h;
        b=(grid.right(stencil(k))-face)/h;
        for power=0:n-1
            M(k,power+1)=face^2*(b^(power+1)-a^(power+1))/(power+1) ...
                +2*face*h*(b^(power+2)-a^(power+2))/(power+2) ...
                +h^2*(b^(power+3)-a^(power+3))/(power+3);
        end
        M(k,:)=M(k,:)/M(k,1);
    end
    unit=zeros(n,1);
    unit(1)=1;
    weights=(M'\unit)';
end

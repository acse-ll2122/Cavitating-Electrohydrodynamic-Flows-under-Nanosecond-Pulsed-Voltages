function grid = mesh1d(C)
    %MESH1D Nonuniform spherical finite-volume grid.
    % Faces follow the prescribed outward-expanding power law. Cell centers are
    % volume centroids, and cell volumes omit the common factor 4*pi.
    N = C.Nr;
    i = (1:N).';
    dr = (C.b - C.a) .* i.^C.meshBeta ./ sum(i.^C.meshBeta);
    grid.faces = [C.a; C.a + cumsum(dr)];
    grid.faces(end) = C.b;
    grid.left = grid.faces(1:end-1);
    grid.right = grid.faces(2:end);
    grid.volume = (grid.right.^3 - grid.left.^3) / 3;
    grid.area = grid.faces.^2;
    grid.centers = (3/4) .* (grid.right.^4 - grid.left.^4) ...
        ./ (grid.right.^3 - grid.left.^3);
    grid.width = grid.right - grid.left;
end

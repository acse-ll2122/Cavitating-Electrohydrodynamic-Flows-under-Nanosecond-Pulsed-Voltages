function op = space_operator(C,grid,state)
    %SPACE_OPERATOR Select B00 central or B01--B04 quadratic/Rusanov fluxes.
    if strcmp(C.spaceMode,'central')
        op=central_operator(C,grid,state);
    else
        op=quadratic_operator(C,grid,state);
    end
end

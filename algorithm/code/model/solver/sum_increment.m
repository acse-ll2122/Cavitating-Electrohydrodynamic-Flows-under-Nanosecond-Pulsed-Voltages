function [value,carry]=sum_increment(value,carry,increment)
    %SUM_INCREMENT TwoSum with a retained low part; no state clipping or filtering.
    % Inputs: high part, retained roundoff, physical increment.
    % Outputs: new high part and its residual low part.
    total=value+increment;
    part=total-value;
    error=(value-(total-part))+(increment-part);
    correction=carry+error;
    value=total+correction;
    carry=correction-(value-total);
end

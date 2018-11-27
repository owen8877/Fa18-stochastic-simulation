function [powIndex, outVal] = pow_out(x, y)
    p = polyfit(x, y, 1);
    powIndex = p(1);
    outVal = p(2);
    fprintf('Pow index seems to be %.2f.\nOutVal: %.2f.\n', powIndex, outVal);
end


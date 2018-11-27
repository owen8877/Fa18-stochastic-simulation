function xi = exp_fit(R, M)
    xi = 1;
    tor = 1e-6;
    maxItr = 1e2;
    step = 1;
    
    fFunc = @(x) matSum((exp(-R/(x^2/4)) - M).^2);
    gFunc = @(x) matSum((exp(-R/(x^2/4)) - M).*exp(-R/(x^2/4)).*R)/x^3;
    
    for itr = 1:maxItr
        g = gFunc(xi); f = fFunc(xi);
%         disp([g f]);
        xi = xi - step * g;
        if f < tor
            break
        end
    end
end

function s = matSum(mat)
    s = sum(sum(mat));
end


function [generator, statLim, distribution, description] = generate_GUE(shallNormalize)
    generator = @generate_GUE_core;
    if shallNormalize
        statLim = 1*3;
        distribution = @(edges) (32/pi^2)*edges.^2.*exp(-edges.^2*4/pi);
    else
        statLim = 4/sqrt(pi)*3;
        distribution = @(edges) (0.5/sqrt(pi))*edges.^2.*exp(-edges.^2/4);
    end
    description = 'GUE';
end

function [H1, H2] = generate_GUE_core(N)
    H1 = generate_1_GUE(N);
    H2 = generate_1_GUE(N);
end

function H = generate_1_GUE(N)
    r = randn(N);
    u = triu(r, 1);
    l = tril(r, -1);
    d = diag(r);
    L = (u + l'*1j)/sqrt(2);
    H = L + L';
    H(1:N+1:N^2) = d;
end
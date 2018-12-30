function [generator, statLim, distribution, description] = generate_Wigner(shallNormalize)
    generator = @generate_Wigner_core;
    if shallNormalize
        statLim = 1*3;
        distribution = @(edges) [];
    else
        statLim = 4/sqrt(pi)*3;
        distribution = @(edges) [];
    end
    description = 'real Wigner';
end

function [H1, H2] = generate_Wigner_core(N)
    r = randn(N);
    u = triu(r, 1);
    l = tril(r, -1);
    H1 = u + u';
    H1(1:N+1:N^2) = randn(N, 1);
    H2 = l + l';
    H2(1:N+1:N^2) = randn(N, 1);
end
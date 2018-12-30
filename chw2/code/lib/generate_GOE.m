function [generator, statLim, distribution, description] = generate_GOE(shallNormalize)
    generator = @generate_GOE_core;
    if shallNormalize
        statLim = 1*3;
        distribution = @(edges) (pi/2)*edges.*exp(-edges.^2*pi/4);
    else
        statLim = sqrt(2*pi)*3;
        distribution = @(edges) (1/4)*edges.*exp(-edges.^2/8);
    end
    description = 'GOE';
end

function [H1, H2] = generate_GOE_core(N)
    r = randn(N);
    u = triu(r, 1);
    l = tril(r, -1);
    % d = diag(r);
    H1 = u + u';
    H1(1:N+1:N^2) = randn(N, 1)*sqrt(2);
    H2 = l + l';
    H2(1:N+1:N^2) = randn(N, 1)*sqrt(2);
end
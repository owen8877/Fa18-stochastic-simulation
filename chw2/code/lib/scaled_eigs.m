function s = scaled_eigs(H, N)
    l = eigs(H, size(H, 1), 0);
    s = sort(l)./sqrt(N);
end
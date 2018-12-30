function s = normalized_eigs(H)
%     l = eigs(@(x) H*x, size(H, 1), size(H, 1), 'smallestreal', 'IsFunctionSymmetric', 1);
    l = eigs(H, size(H, 1), 0);
    s = sort(l);
    s = s(2:end)-s(1:end-1);
end
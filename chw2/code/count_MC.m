function crossResults = count_MC(edges, crossTimes, allTries, N, generator)
    crossResults = zeros(crossTimes, numel(edges)-1);
    for j = 1:crossTimes
        results = zeros(allTries, N);
        for tries = 2:2:allTries
            [H1, H2] = generator(N);
            results(tries-1, :) = scaled_eigs(H1, N)';
            results(tries, :) = scaled_eigs(H2, N)';
        end
        
        crossResults(j, :) = histcounts(reshape(results, [numel(results) 1]), edges);
        fprintf('Finished round %d.\n', j)
    end
end
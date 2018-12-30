function crossResults = main_MC(edges, crossTimes, allTries, N, generator, shallNormalize)
    crossResults = zeros(crossTimes, numel(edges)-1);
    for j = 1:crossTimes
        results = zeros(allTries, N-1);
        for tries = 2:2:allTries
            [H1, H2] = generator(N);
            results(tries-1, :) = normalized_eigs(H1)';
            results(tries, :) = normalized_eigs(H2)';
        end

        if shallNormalize
            results = results ./ mean(results, 'all');
        end
        crossResults(j, :) = histcounts(reshape(results, [numel(results) 1]), edges);
        fprintf('Finished round %d.\n', j)
    end
end
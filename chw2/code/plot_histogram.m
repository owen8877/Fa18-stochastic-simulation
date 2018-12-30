clear; %clc
addpath lib

allTries = 1e4; crossTimes = 10;
Mbins = 50;
N = 10;
descriptions = {'GOE', 'GUE', 'real Wigner'};

%% loading
results = cell(numel(descriptions), 1);
for k = 1:numel(descriptions)
    results{k} = load(sprintf('data/N%d-allT%d-crossT%d-hist-%s', ...
        N, allTries, crossTimes, descriptions{k}));
end

%% For GOE; before-after normalization
figure(5); hold on
xlim([-2.5 2.5]), ylim([0 0.35]); grid on
helper(results{1}, 'GOE')
helper(results{2}, 'GUE')
helper(results{3}, 'real Wigner')
edges = results{1}.edges;
plot(edges, sqrt(4-edges.^2)/2/pi, '-.', 'DisplayName', 'theory');
legend; title('Distribution of eigenvalues (N=10)')
xlabel('\lambda/\sqrt{N}', 'interpreter', 'latex'); ylabel('frequency')
set(gcf, 'Position', [488 540 630 252])

function helper(data, tag)
    allTries = data.allTries;
    crossResults = data.crossResults;
    edges = data.edges;
    N = data.N;
    scaling = 1/allTries/(edges(2)-edges(1))/N;
    errorbar((edges(2:end)+edges(1:end-1))/2, mean(crossResults, 1)*scaling, 2*std(crossResults, 1)*scaling, 'DisplayName', tag);
end
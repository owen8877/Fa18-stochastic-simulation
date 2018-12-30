clear; %clc
addpath lib

N = 2;
allTries = 1e4; crossTimes = 10;
Mbins = 50;
shallNormalize = true;
g = @generate_GOE;

[generator, statLim, distribution, description] = g(shallNormalize);
edges = (0:1/Mbins:1)*statLim;

tic
crossResults = main_MC(edges, crossTimes, allTries, N, generator, shallNormalize);
toc

if N == 2
    desired_d = distribution(edges);
else
    desired_d = [];
end

save(sprintf('data/N%d-allT%d-crossT%d-norm%d-%s', ...
    N, allTries, crossTimes, shallNormalize, description), ...
    'N', 'allTries', 'crossTimes', 'edges', 'crossResults', 'desired_d')
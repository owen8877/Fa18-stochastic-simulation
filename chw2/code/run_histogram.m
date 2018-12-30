clear; %clc
addpath lib

N = 10;
allTries = 1e4; crossTimes = 10;
Mbins = 50;
g = @generate_GOE;

[generator, ~, ~, description] = g(true);
edges = -2.5:5/Mbins:2.5;

tic
crossResults = count_MC(edges, crossTimes, allTries, N, generator);
toc

save(sprintf('data/N%d-allT%d-crossT%d-hist-%s', ...
    N, allTries, crossTimes, description), ...
    'N', 'allTries', 'crossTimes', 'edges', 'crossResults')
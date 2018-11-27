clear; %clc
addpath lib

folder = 'q54';
check_folder(['./output/' folder]);

% query.N = 12;
% query.sample = 800;
% query.preWait = 200;
% query.outputInt = 10000;

query.N = 14;
query.sample = 1200;
query.preWait = 400;
query.outputInt = 20000;
query.J = 1.0;
query.h = 0.0;
query.dimension = 3;
query.allTries = 5;

load(sprintf('data/beta-all-%dd.mat', query.dimension), 'betas')

for beta = betas
    query_run = query;
    query_run.beta = beta;
    runTest(folder, query_run);
end

c = clock();
check_folder(['./info/' folder]);
save(sprintf('./info/%s/%d-betax%d-%d-%d.%d.mat', ...
    folder, query.N, numel(betas), c(3), c(4), c(5)), ...
    'query', 'betas', 'times')
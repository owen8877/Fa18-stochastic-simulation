clear; %clc
addpath lib

folder = 'q2';
check_folder(['./output/' folder]);

query.N = 30; query.sample = 1600; query.preWait = 600; query.outputInt = 10000;
query.J = 1.0;
query.dimension = 2;
query.allTries = 10;
query.noDump = true;

load(sprintf('data/beta-fine-%dd.mat', query.dimension), 'betas')
% load(sprintf('data/beta-all-%dd.mat', query.dimension), 'betas')

hs = 0.000;
betas = [0.36 0.39 0.42 0.45 0.48 0.51];
times = [];
for beta = betas
    for h = hs
        query_run = query;
        query_run.beta = beta;
        query_run.h = h;

        query_run.preWait = round(beta^(1/4) * query.preWait);
        query_run.sample = round(beta^(1/4) * query.sample);
        t = runTest(folder, query_run);
        times = [times; t];
    end
end

c = clock();
check_folder(['./info/' folder]);
save(sprintf('./info/%s/%d-betax%d-%d-%d.%d.mat', ...
    folder, query.N, numel(betas), c(3), c(4), c(5)), ...
    'query', 'betas', 'times', 'hs')
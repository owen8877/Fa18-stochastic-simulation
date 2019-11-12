clear; %clc
addpath lib

folder = 'q7';
check_folder(['./output/' folder]);

% query.N = 15; query.sample = 750; query.preWait = 150; query.outputInt = 4500;
% query.N = 20; query.sample = 1000; query.preWait = 200; query.outputInt = 6000;
query.N = 25; query.sample = 1250; query.preWait = 250; query.outputInt = 6000;
query.J = 1.0;
query.h = 0.0;
query.dimension = 2;
query.allTries = 5;

load(sprintf('data/beta-all-%dd.mat', query.dimension), 'betas')

times = [];
for beta = betas
    query_run = query;
    query_run.beta = beta;
    
    query_run.preWait = round(beta^(1/4) * query.preWait);
    query_run.sample = round(beta^(1/4) * query.sample);
    t = runTestbcc(folder, query_run);
    times = [times; t];
end

c = clock();
check_folder(['./info/' folder]);
save(sprintf('./info/%s/%d-betax%d-%d-%d.%d.mat', ...
    folder, query.N, numel(betas), c(3), c(4), c(5)), ...
    'query', 'betas', 'times')
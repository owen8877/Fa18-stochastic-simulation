clear; %clc
addpath lib

folder = 'q1';
check_folder(['./output/' folder]);

% query.N = 20; query.sample = 800; query.preWait = 200; query.outputInt = 5000;
% query.N = 25; query.sample = 1200; query.preWait = 400; query.outputInt = 10000;
query.N = 30; query.sample = 1600; query.preWait = 600; query.outputInt = 10000;
% query.N = 35; query.sample = 2400; query.preWait = 800; query.outputInt = 10000;
% query.N = 40; query.sample = 2400; query.preWait = 800; query.outputInt = 15000;
query.J = 1.0;
query.h = 0.0;
query.dimension = 2;
query.allTries = 5;

% load(sprintf('data/beta-fine-%dd.mat', query.dimension), 'betas')
load(sprintf('data/beta-all-%dd.mat', query.dimension), 'betas')

times = [];
for beta = betas
    query_run = query;
    query_run.beta = beta;
    
    query_run.preWait = round(beta^(1/4) * query.preWait);
    query_run.sample = round(beta^(1/4) * query.sample);
%     query_run.preWait = round(query.preWait);
%     query_run.sample = round(query.sample);
    t = runTest(folder, query_run);
    times = [times; t];
end

c = clock();
check_folder(['./info/' folder]);
save(sprintf('./info/%s/%d-betax%d-%d-%d.%d.mat', ...
    folder, query.N, numel(betas), c(3), c(4), c(5)), ...
    'query', 'betas', 'times')
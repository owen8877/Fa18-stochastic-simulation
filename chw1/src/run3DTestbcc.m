clear; %clc
addpath lib

folder = 'q8bcc';
check_folder(['./output/' folder]);

configs = [ ...
     8,  600,  120,  8000; ...
    10,  600,  150, 10000; ...
%   12,  800,  200, 10000; ...
    14, 1200,  200, 10000; ...
    16, 1200,  200, 12500; ...
    18, 1200,  200, 15000; ...
];

% query.N = 12; query.sample = 800;
% query.preWait = 200; query.outputInt = 10000;

query.J = 1.0;
query.h = 0.0;
query.dimension = 3;
query.allTries = 5;

load(sprintf('data/beta-fine-%dd-bcc.mat', query.dimension), 'betas')

for co = configs'
    query.N = co(1);
    query.sample = co(2);
    query.preWait = co(3);
    query.outputInt = co(4);
    for beta = betas
        query_run = query;
        query_run.beta = beta;
        times = runTestbcc(folder, query_run);
    end

    c = clock();
    check_folder(['./info/' folder]);
    save(sprintf('./info/%s/%d-betax%d-%d-%d.%d.mat', ...
        folder, query.N, numel(betas), c(3), c(4), c(5)), ...
        'query', 'betas', 'times')
end
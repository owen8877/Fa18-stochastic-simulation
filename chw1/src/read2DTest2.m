clear; %clc
addpath lib

folder = 'q2';

query.N = 30;
query.J = 1.0;
hs = 0.000:0.005:0.030;
query.dimension = 2;
query.allTries = 10;

betas = [0.36 0.39 0.42 0.45 0.48 0.51];
betaN = numel(betas);

hf2 = figure(2);
ax2o = axes();
ax2 = axes('parent', hf2, 'position', [0.45 0.35 0.4 0.35]);

[magMeans, magVars] = readTestAlongHs(folder, query, betas, hs);

figure(2); clf; hold on
for i = 1:numel(betas)
    beta = betas(i);
    blockData = abs(magMeans(:, :, i))/query.N^2;
    errorbar(hs, mean(blockData, 1), 2*std(blockData, 1), 's-', 'DisplayName', num2str(beta))
end
title(legend(), '\beta')
ylim([0 1]); grid on
xlabel('External field intensity h')
ylabel('m')

%% Helpers
function [magMeans, magVars] = readTestAlongHs(qn, query, betas, hs)
    addpath lib
    
    N = query.N;
    dimension = query.dimension;
    allTries = query.allTries;

    hsN = numel(hs);
    betasN = numel(betas);
    magMeans = nan*ones(allTries, hsN, betasN);
    magVars = nan*ones(allTries, hsN, betasN);
    
    for i = 1:betasN
        beta = betas(i);
        for l = 1:hsN
            h = hs(l);
            for j = 1:allTries
                saveFileName = sprintf('./output/%s/%04.4f-%.4f-%d-%d.mat.%d', ...
                    qn, beta, h, N, dimension, j);
                if ~exist(saveFileName, 'file')
                    continue
                end
                data = load(saveFileName, '-ascii');

                data_mag = data(:, 2);
                magMeans(j, l, i) = mean(data_mag);
                magVars(j, l, i) = var(data_mag);
            end
        end
    end
end
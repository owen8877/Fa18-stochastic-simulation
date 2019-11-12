clear; %clc
addpath lib

folder = 'q8bcc';

Ns = [12];
query.J = 1.0;
query.h = 0.0;
query.sample = 800;
query.preWait = 200;
% query.outputInt = 5000;
query.dimension = 3;
query.allTries = 5;

load(sprintf('data/beta-fine-%dd-bcc.mat', query.dimension), 'betas')
betaN = numel(betas);
% T_critical = 2.27;
% idx = 15;
T_critical = 6;
idx = 7;

for N = Ns
    query.N = N;
    [means, vars, magMeans, magVars, xiMeans, corrMat] = readTestbcc(folder, query, betas);

    figure(1)
    subplot(2, 2, 1); part1plotHelper(betas, means, query, true, 'Inner Energy', num2str(query.N))
    subplot(2, 2, 2); part1plotHelper(betas, betas.^2 .* vars, query, true, 'Specific Heat', num2str(query.N))
%     subplot(2, 3, 3); part1plotHelper(betas, xiMeans, query, false, 'Characterisitic Length', num2str(query.N))
    subplot(2, 2, 3); part1plotHelper(betas, magMeans, query, true, 'Magneti...', num2str(query.N))
    subplot(2, 2, 4); part1plotHelper(betas, betas.^2 .* magVars, query, true, 'Magneti Kai...', num2str(query.N))
    
%     if N == Ns(end)
%         % correlation length varies with distance and beta
%         figure(3);
%         [disX, disY] = meshgrid(1:N, betas);
%         mesh(disX, disY, corrMat');
%         xlabel('r');
%         ylabel('beta')
%     end

    figure(2)
    subplot(2, 2, 1); part2plotHelper(@loglog, idx, 'high', betas, betas.^2 .* vars, T_critical, query, true, 'Specific Heat', num2str(query.N))
    subplot(2, 2, 2); part2plotHelper(@semilogx, idx, 'low', betas, betas.^2 .* vars, T_critical, query, true, 'Specific Heat', num2str(query.N))
    subplot(2, 2, 3); part2plotHelper(@loglog, idx, 'low', betas, magMeans, T_critical, query, true, 'Magneti...', num2str(query.N))
%     subplot(2, 2, 4); part2plotHelper(@loglog, idx, 'high', betas, xiMeans, T_critical, query, false, 'Characterisitic Length', num2str(query.N))
end

function part1plotHelper(betas, data, query, dallN, title_txt, label_txt)
    if dallN    
        allN = query.N^query.dimension;
    else
        allN = 1;
    end
    hold on
    errorbar(betas, ...
        mean(data, 1) / allN, ...
        2*std(data, 1) ./ (sqrt(query.sample*query.allTries) * allN), ...
        'DisplayName', label_txt);
    title(title_txt)
    xlabel('\beta')
    grid on
    legend
end

function part2plotHelper(plotFunc, idx, rangeType, betas, data, T_critical, query, dallN, title_txt, label_txt)
    signD = 1;
    switch rangeType
        case 'low'
            betas = betas(idx+2:end);
            data = data(:, idx+2:end);
            title_app = '; T<T*';
        case 'high'
            betas = betas(1:idx-1);
            data = data(:, 1:idx-1);
            signD = -1;
            title_app = '; T>T*';
        otherwise
            title_app = '; all range T';
    end
    if dallN    
        allN = query.N^query.dimension;
    else
        allN = 1;
    end
    temp = 1./betas;
    xData = (1-temp./T_critical)*signD;
    yData = mean(data, 1)/allN;
    eData = std(data, 1) ./ (sqrt(query.sample*query.allTries) * allN);
    hold on
	errorbar(xData, yData, eData, 'DisplayName', label_txt);
    title([title_txt title_app]);
    grid on
    legend
    
    if strcmp(func2str(plotFunc), 'loglog')
        p = polyfit(log(xData), log(yData), 1);
        set(gca, 'XScale', 'log')
        set(gca, 'YScale', 'log')
    elseif strcmp(func2str(plotFunc), 'semilogx')
        p = polyfit(log(xData), yData, 1);
        set(gca, 'XScale', 'log')
        set(gca, 'YScale', 'linear')
    elseif strcmp(func2str(plotFunc), 'semilogy')
        p = polyfit(xData, log(yData), 1);
        set(gca, 'XScale', 'linear')
        set(gca, 'YScale', 'log')
    end
    
    fprintf('Slope for %s%s: %.4f\n', title_txt, title_app, p(1));
end
clear; %clc
addpath lib

folder = 'q54';

Ns = [8 10 12 14];
query.J = 1.0;
query.h = 0.0;
query.dimension = 3;
query.allTries = 5;

load(sprintf('data/beta-all-%dd.mat', query.dimension), 'betas')
betaN = numel(betas);
T_critical = 4.425;
idx = 12;

hf2 = figure(2);
ax2o = axes();
ax2 = axes('parent', hf2, 'position', [0.25 0.45 0.2 0.35]);

for N = Ns
    if N < 13
        continue
    end
    query.N = N; fprintf('\nNow N=%d\n', N)
    [means, vars, magMeans, magVars, xiMeans, corrMat] = readTest(folder, query, betas);
    
    figure(1)
    part1plotHelper(betas, means, query, true, 'Inner Energy', num2str(query.N))
    xlabel('Temperature T')
    ylabel('Internal Energy u')
    
    figure(2);
    axes(ax2o);
    part1plotHelper(betas, betas.^2 .* vars, query, true, 'Specific Heat', num2str(query.N))
    xlabel('Temperature T')
    ylabel('Specific heat c')
    axes(ax2);
    box on
    part1plotHelper(betas, betas.^2 .* vars, query, true, '', num2str(query.N))
    xlim([4.2 4.55]); ylim([1.4 2.5])
    
    figure(3)
    subplot(2, 2, 1)
    part2plotHelper(@loglog, idx, 'high', betas, betas.^2 .* vars, T_critical, query, true, 'Specific Heat', num2str(query.N))
    xlabel('\epsilon')
    ylabel('Specific heat c (Log scale)')
    subplot(2, 2, 2)
    part2plotHelper(@semilogx, idx, 'high', betas, betas.^2 .* vars, T_critical, query, true, 'Specific Heat', num2str(query.N))
    xlabel('\epsilon')
    ylabel('Specific heat c (Semilog-x scale)')
    subplot(2, 2, 3)
    part2plotHelper(@loglog, idx, 'low', betas, betas.^2 .* vars, T_critical, query, true, 'Specific Heat', num2str(query.N))
    xlabel('\epsilon')
    ylabel('Specific heat c (Log scale)')
    subplot(2, 2, 4)
    part2plotHelper(@semilogx, idx, 'low', betas, betas.^2 .* vars, T_critical, query, true, 'Specific Heat', num2str(query.N))
    xlabel('\epsilon')
    ylabel('Specific heat c (Semilog-x scale)')
    
    figure(4)
    subplot(1, 2, 1)
    part1plotHelper(betas, xiMeans, query, false, 'Characterisitic Length', num2str(query.N))
    set(gca, 'YScale', 'log')
    xlabel('T')
    ylabel('\xi')
    subplot(1, 2, 2)
    part2plotHelper(@loglog, idx, 'high', betas, xiMeans, T_critical, query, false, 'Characterisitic Length', num2str(query.N))
    xlabel('\epsilon')
    ylabel('\xi (Log scale)')
    
    figure(5)
    subplot(1, 3, 1)
    part1plotHelper(betas, magMeans, query, true, 'Magnetization', num2str(query.N))
    xlabel('T'); xlim([1 6.25])
    ylabel('m')
    subplot(1, 3, 2)
    part2plotHelper(@loglog, idx, 'low', betas, magMeans, T_critical, query, true, 'Magnetization', num2str(query.N))
    xlabel('\epsilon')
    ylabel('m (Log Scale)')
    subplot(1, 3, 3)
    part2plotHelperForM(@semilogy, idx, 'low', betas, magMeans, T_critical, query, true, 'Magnetization', num2str(query.N))
    xlabel('\beta')
    ylabel('m (Semilog-y Scale)')
    
    figure(6)
    subplot(1, 2, 1)
    part1plotHelper(betas, betas .* magVars, query, true, 'Magnetic Susceptibility', num2str(query.N))
    xlabel('T')
    ylabel('\chi')
    subplot(1, 2, 2)
    part2plotHelper(@loglog, idx, 'high', betas, betas .* magVars, T_critical, query, true, 'Magnetic Susceptibility', num2str(query.N))
    xlabel('\epsilon')
    ylabel('\chi (Log Scale)')
end
figure(1); hold on;                               plot([1 1]*T_critical, ylim, 'k:', 'DisplayName', 'Critical Temperature')
figure(2); axes(ax2o); hold on;                   plot([1 1]*T_critical, ylim, 'k:', 'DisplayName', 'Critical Temperature')
figure(2); axes(ax2); hold on;                    plot([1 1]*T_critical, ylim, 'k:', 'DisplayName', 'Critical Temperature')
figure(4); s1 = subplot(1, 2, 1); hold(s1, 'on'); plot([1 1]*T_critical, ylim, 'k:', 'DisplayName', 'Critical Temperature')
figure(5); s1 = subplot(1, 3, 1); hold(s1, 'on'); plot([1 1]*T_critical, ylim, 'k:', 'DisplayName', 'Critical Temperature')
figure(6); s1 = subplot(1, 2, 1); hold(s1, 'on'); plot([1 1]*T_critical, ylim, 'k:', 'DisplayName', 'Critical Temperature')

figure(3); s1 = subplot(2, 2, 1); hold(s1, 'on');
SIHelper(@loglog, idx, 'high', betas, T_critical, -0.8803, -3.0039)
figure(3); s2 = subplot(2, 2, 2); hold(s2, 'on');
SIHelper(@semilogx, idx, 'high', betas, T_critical, -0.4595, -0.5259)
figure(3); s3 = subplot(2, 2, 4); hold(s3, 'on');
SIHelper(@semilogx, idx, 'low', betas, T_critical, -0.5939, -0.2074)
figure(4); s2 = subplot(1, 2, 2); hold(s2, 'on');
SIHelper(@loglog, idx, 'high', betas, T_critical, -0.2264, -0.7334)
figure(5); s2 = subplot(1, 3, 2); hold(s2, 'on');
SIHelper(@loglog, idx, 'high', betas, T_critical, 0.2006, 0.1759); ylim([0.6 1])
figure(5); s3 = subplot(1, 3, 3); hold(s3, 'on');
plot(betas(idx+1:end), exp(-11.4506 * betas(idx+1:end) + 0.3263), 'd-.', 'DisplayName', 'fit'); ylim([1e-5 4e-1])
figure(6); s2 = subplot(1, 2, 2); hold(s2, 'on');
SIHelper(@loglog, idx, 'high', betas, T_critical, -1.2070, -1.1522)

positions = [ ...
    488.2000 532.2000 555.2000 252.0000; ...
    488.2000 516.2000 715.2000 268.0000; ...
    137.0000 101.0000 731.2000 524     ; ...
    94.6000  321.8000 907.2000 340.0000; ...
    138.6    310.6    1011.2   292.0   ; ...
    94.6000  321.8000 907.2000 340.0000; ...
];
for i = 1:6
    figure(i); set(gcf, 'Position', positions(i, :))
end

%% Peak Points Approximation
% peakPoints = [];
% peakVals = [];
% for N = Ns
%     query.N = N;
%     [means, vars, magMeans, magVars, xiMeans, corrMat] = readTest(folder, query, betas);
%     
%     cs = mean(betas.^2 .* vars, 1)/N^2;
%     p = polyfit(1./betas(10:18), cs(10:18), 2);
%     peakPoints = [peakPoints p(2)/(-2*p(1))];
%     peakVals = [peakVals p(3)-p(2)^2/(4*p(1))];
% end
% [pp1, pp2] = pow_out(1./Ns, peakPoints);
% [pv1, pv2] = pow_out(1./Ns, peakVals);
% figure(7)
% s1 = subplot(1, 2, 1); hold(s1, 'on')
% plot(1./Ns, peakPoints, '-s', 'DisplayName', 'Exp');
% plot([0 1/Ns(1)], [0 1/Ns(1)]*pp1+pp2, '--', 'DisplayName', 'Fit')
% title('Fit for peak point.')
% legend; grid on
% xlabel('1/N'); ylabel('T_c(N)')
% s2 = subplot(1, 2, 2); hold(s2, 'on')
% plot(1./Ns, peakVals, '-s', 'DisplayName', 'Exp');
% plot([0 1/Ns(1)], [0 1/Ns(1)]*pv1+pv2, '--', 'DisplayName', 'Fit')
% title('Fit for peak value.')
% legend; grid on
% xlabel('1/N'); ylabel('c_{max}(N)')

%% Helpers
function SIHelper(plotFunc, idx, rangeType, betas, T_critical, slope, intercept)
    signD = 1;
    switch rangeType
        case 'low'
            betas = betas(idx+2:end);
        case 'high'
            betas = betas(1:idx-1);
            signD = -1;
    end
    temp = 1./betas;
    xData = (1-temp./T_critical)*signD;
    switch func2str(plotFunc)
        case 'semilogx'
            yData = slope * log(xData) + intercept;
        case 'semilogy'
            yData = exp(slope * xData + intercept);
        case 'loglog'
            yData = exp(slope * log(xData) + intercept);
    end
    plot(xData, yData, 'd-.', 'DisplayName', 'fit')
end

function part1plotHelper(betas, data, query, dallN, title_txt, label_txt)
    if dallN    
        allN = query.N^query.dimension;
    else
        allN = 1;
    end
    hold on
    errorbar(1./betas, ...
        mean(data, 1) / allN, ...
        2*std(data, 1) ./ (sqrt(size(data, 2)*query.allTries) * allN), ...
        'DisplayName', label_txt);
    title(title_txt)
    grid on
    title(legend, 'N')
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
    eData = std(data, 1) ./ (sqrt(size(data, 2)*query.allTries) * allN);
    hold on
	errorbar(xData, yData, eData, 'DisplayName', label_txt);
    title([title_txt title_app]);
    grid on
    title(legend, 'N')
    
    nanid = isnan(yData);
    xData_ = xData(~nanid);
    yData_ = yData(~nanid);
    if strcmp(func2str(plotFunc), 'loglog')
        p = polyfit(log(xData_), log(yData_), 1);
        set(gca, 'XScale', 'log')
        set(gca, 'YScale', 'log')
    elseif strcmp(func2str(plotFunc), 'semilogx')
        p = polyfit(log(xData_), yData_, 1);
        set(gca, 'XScale', 'log')
        set(gca, 'YScale', 'linear')
    elseif strcmp(func2str(plotFunc), 'semilogy')
        p = polyfit(xData_, log(yData_), 1);
        set(gca, 'XScale', 'linear')
        set(gca, 'YScale', 'log')
    end
    
    fprintf('%s%s: slope=%.4f, intercept=%.4f\n', title_txt, title_app, p(1), p(2));
end

function part2plotHelperForM(plotFunc, idx, ~, betas, data, ~, query, dallN, title_txt, label_txt)
    betas = betas(idx+2:end);
    data = data(:, idx+2:end);
    title_app = '; T<T*';
    
    if dallN    
        allN = query.N^query.dimension;
    else
        allN = 1;
    end
    xData = betas;
    yData = 1-mean(data, 1)/allN;
    eData = std(data, 1) ./ (sqrt(size(data, 2)*query.allTries) * allN);
    hold on
	errorbar(xData, yData, eData, 'DisplayName', label_txt);
    title([title_txt title_app]);
    grid on
    title(legend, 'N')
    
    nanid = isnan(yData);
    xData_ = xData(~nanid);
    yData_ = yData(~nanid);
    if strcmp(func2str(plotFunc), 'loglog')
        p = polyfit(log(xData_), log(yData_), 1);
        set(gca, 'XScale', 'log')
        set(gca, 'YScale', 'log')
    elseif strcmp(func2str(plotFunc), 'semilogx')
        p = polyfit(log(xData_), yData_, 1);
        set(gca, 'XScale', 'log')
        set(gca, 'YScale', 'linear')
    elseif strcmp(func2str(plotFunc), 'semilogy')
        p = polyfit(xData_, log(yData_), 1);
        set(gca, 'XScale', 'linear')
        set(gca, 'YScale', 'log')
    end
    
    fprintf('%s%s: slope=%.4f, intercept=%.4f\n', title_txt, title_app, p(1), p(2));
end
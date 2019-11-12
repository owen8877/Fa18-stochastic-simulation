function [means, vars, magMeans, magVars, xiMeans, corrMat] = readTestbcc(qn, query, betas)
    addpath lib
    
    N = query.N;
    dimension = query.dimension;
    h = query.h;
    allTries = query.allTries;

    betaN = numel(betas);
    
    processedName = sprintf('./processed/%s/%d-%d-%d.mat.%d', ...
        qn, betaN, N, dimension, allTries);
    if exist(processedName, 'file')
        % load(processedName, '-mat', 'means', 'vars', 'magMeans', 'magVars', 'xiMeans', 'corrMat')
        load(processedName, '-mat', 'means', 'vars', 'magMeans', 'magVars')
        xiMeans = nan*ones(allTries, betaN);
        corrMat = nan*ones(N, betaN);
        return
    end
            
    means = nan*ones(allTries, betaN);
    vars = nan*ones(allTries, betaN);
    magMeans = nan*ones(allTries, betaN);
    magVars = nan*ones(allTries, betaN);
    xiMeans = nan*ones(allTries, betaN);
    corrMat = nan*ones(N, betaN);

%     if mod(N, 2) == 0
%         d = disMatG(N/2+1, N/2+1);
%         rs = [d fliplr(d(:, 2:end-1))];
%         rs = [rs; flipud(rs(2:end-1, :))];
%     else
%         d = disMatG(N/2+1/2, N/2+1/2);
%         rs = [d fliplr(d(:, 2:end))];
%         rs = [rs; flipud(rs(2:end, :))];
%     end
    
    for i = 1:betaN
        %　corrTmp = zeros(N, N, allTries);
        for j = 1:allTries
            beta = betas(i);
            saveFileName = sprintf('./output/%s/%04.4f-%.4f-%d-%d.mat.%d', ...
                qn, beta, h, N, dimension, j);
            if ~exist(saveFileName, 'file')
                continue
            end
            data = load(saveFileName, '-ascii');
            
            data_energy = data(:, 1);
            data_mag = data(:, 2);
            % acfArr = acf(data_energy, 10);

            % corrDump = sprintf('./output/%s/%04.4f-%.4f-%d-%d.dump.%d', ...
            %     qn, beta, h, N, dimension, j);
            % corrData = load(corrDump, '-ascii');
            % corrTmp(:, :, j) = corrData;

            means(j, i) = mean(data_energy);
            vars(j ,i) = var(data_energy);
            magMeans(j, i) = mean(abs(data_mag));
            magVars(j, i) = var(data_mag);
            % xiMeans(j, i) = sum(sum(corrData.*rs))/(sample*N^dimension/2/sqrt(pi));
            % xiMeans(j, i) = exp_fit(rs, corrData);
        end
        % mcT = mean(corrTmp, 3);
        % xiMeans(:, i) = exp_fit(rs, mcT);
        % ydata = (mcT(:, 1) + mcT(1, :)')/2/size(data, 1);
        % xdata = min(0:N-1, N:-1:1)';
        % corrMat(:, i) = ydata;
        % xiMeans(:, i) = LMFnlsq(@(x) exp(-xdata/x)-ydata, 1);
        
        % correlation figure
        % if mod(i, 6) == 0 && beta < 0.8 && N == 14
        %     figure(20); hold on
        %     plot(xdata(1:floor(end/2)), ydata(1:floor(end/2)), 's-', 'DisplayName', num2str(beta))
        % end

        % figure(10); title(['beta=' num2str(beta)])
        % subplot(2, 2, 1); plot(1:10, acfArr); grid on
        % subplot(2, 2, 2);
        % mesh(circshift(mean(corrTmp, 3)/sample, [N/2, N/2]));
        % title(num2str(cesseMeans(end, i)));
        % subplot(2, 2, 3); plot(data(:, 1))
        % subplot(2, 2, 4); plot(data(:, 2))
        % pause
    end
    
    % save(processedName, 'means', 'vars', 'magMeans', 'magVars', 'xiMeans', 'corrMat')
    save(processedName, 'means', 'vars', 'magMeans', 'magVars')
end

function disMat = disMatG(r, c)
    disMat = sqrt((0:c-1).^2 + ((0:r-1)').^2);
end
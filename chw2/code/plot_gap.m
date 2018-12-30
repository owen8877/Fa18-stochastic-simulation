clear; %clc
addpath lib

allTries = 1e4; crossTimes = 10;
Mbins = 50;
Ns = [2 4 10];
descriptions = {'GOE', 'GUE', 'real Wigner'};

%% loading
results = cell(numel(Ns), 2, numel(descriptions));
for i = 1:numel(Ns)
    for j = 1:2
        for k = 1:numel(descriptions)
            results{i, j, k} = load(sprintf('data/N%d-allT%d-crossT%d-norm%d-%s', ...
                Ns(i), allTries, crossTimes, j-1, descriptions{k}));
        end
    end
end

%% For GOE; before-after normalization
figure(1)
s1 = subplot(1, 2, 1); hold(s1, 'on'); ylim([0 0.7]); xlim([0 6]); grid on
helper(results{1, 1, 1}, '2', true, '2 (theory)')
helper(results{2, 1, 1}, '4', false, '')
helper(results{3, 1, 1}, '10', false, '')
legend; title('Gap of eigenvalues (GOE, original)')
xlabel('|\lambda_1-\lambda_2|'); ylabel('frequency')

s2 = subplot(1, 2, 2); hold(s2, 'on'); ylim([0 0.9]); xlim([0 3]); grid on
helper(results{1, 2, 1}, '2', true, '2 (theory)')
helper(results{2, 2, 1}, '4', false, '')
helper(results{3, 2, 1}, '10', false, '')
legend; title('Gap of eigenvalues (GOE, normalized)')
xlabel('|\lambda_1-\lambda_2|/\langle|\lambda_1-\lambda_2|\rangle'); ylabel('frequency')

set(gcf, 'Position', [301 313 851 300])

%% For GUE; before-after normalization
figure(2)
s1 = subplot(1, 2, 1); hold(s1, 'on'); ylim([0 0.9]); xlim([0 6]); grid on
helper(results{1, 1, 2}, '2', true, '2 (theory)')
helper(results{2, 1, 2}, '4', false, '')
helper(results{3, 1, 2}, '10', false, '')
legend; title('Gap of eigenvalues (GUE, original)')
xlabel('|\lambda_1-\lambda_2|'); ylabel('frequency')

s2 = subplot(1, 2, 2); hold(s2, 'on'); ylim([0 1.1]); xlim([0 3]); grid on
helper(results{1, 2, 2}, '2', true, '2 (theory)')
helper(results{2, 2, 2}, '4', false, '')
helper(results{3, 2, 2}, '10', false, '')
legend; title('Gap of eigenvalues (GUE, normalized)')
xlabel('|\lambda_1-\lambda_2|/\langle|\lambda_1-\lambda_2|\rangle'); ylabel('frequency')

set(gcf, 'Position', [301 313 851 300])

%% For real Wigner; before-after normalization
figure(3)
s1 = subplot(1, 2, 1); hold(s1, 'on'); ylim([0 0.7]); xlim([0 6]); grid on
helper(results{1, 1, 3}, '2', false, '')
helper(results{2, 1, 3}, '4', false, '')
helper(results{3, 1, 3}, '10', false, '')
legend; title('Gap of eigenvalues (real Wigner, original)')
xlabel('|\lambda_1-\lambda_2|'); ylabel('frequency')

s2 = subplot(1, 2, 2); hold(s2, 'on'); ylim([0 0.9]); xlim([0 3]); grid on
helper(results{1, 2, 3}, '2', false, '')
helper(results{2, 2, 3}, '4', false, '')
helper(results{3, 2, 3}, '10', false, '')
legend; title('Gap of eigenvalues (real Wigner, normalized)')
xlabel('|\lambda_1-\lambda_2|/\langle|\lambda_1-\lambda_2|\rangle'); ylabel('frequency')

set(gcf, 'Position', [301 313 851 300])

%% three comparison
figure(4); hold on
ylim([0 1.1]); xlim([0 3]); grid on
helper(results{3, 2, 1}, 'GOE', false, '')
helper(results{3, 2, 2}, 'GUE', false, '')
helper(results{3, 2, 3}, 'real Wigner', false, '')
legend; title('Gap of eigenvalues (N=10, normalized)')
xlabel('|\lambda_1-\lambda_2|/\langle|\lambda_1-\lambda_2|\rangle'); ylabel('frequency')
set(gcf, 'Position', [488 540 595 244])

function helper(data, tag, theoryOn, tag2)
    allTries = data.allTries;
    crossResults = data.crossResults;
    edges = data.edges;
    N = data.N;
    scaling = 1/allTries/(edges(2)-edges(1))/(N-1);
    errorbar((edges(2:end)+edges(1:end-1))/2, mean(crossResults, 1)*scaling, 2*std(crossResults, 1)*scaling, 'DisplayName', tag);
    if theoryOn && numel(data.desired_d) > 0
        plot(edges, data.desired_d, '-.', 'DisplayName', tag2);
    end
end
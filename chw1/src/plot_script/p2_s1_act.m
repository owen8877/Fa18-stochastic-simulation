clear; %clc

acf1 = acfHelper(0.48);
acf2 = acfHelper(0.435);
acf3 = acfHelper(0.4);

figure
semilogy(0:10, [1 acf1'], 'DisplayName', '0.48'); hold on
semilogy(0:10, [1 acf2'], 'DisplayName', '0.435');
semilogy(0:10, [1 acf3'], 'DisplayName', '0.4');
xlabel('jump interval')
ylabel('<H(\sigma_t)H(\sigma_{t+\Deltat})>')
title(legend, '\beta'); grid on

function acf1 = acfHelper(beta)
    qn = 'q1';
    h = 0.0;
    N = 30;
    dimension = 2;
    j = 1;

    saveFileName = sprintf('../output/%s/%04.4f-%.4f-%d-%d.mat.%d', ...
                    qn, beta, h, N, dimension, j);

    data = load(saveFileName, '-ascii');
    data_energy = data(:, 1);
    data_mag = data(:, 2);

    acf1 = abs(acf(data_energy, 10));
end
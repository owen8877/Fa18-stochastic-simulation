clear; %clc

q.N = 50;
q.sample = 1;
q.preWait = 1000;
q.outputInt = 1e4;
q.dimension = 2;

q.J = 1.0;
q.h = 0.0;

betas = [0.2 0.4407 0.7];

figure(15);
for i = 1:numel(betas)
    q.beta = betas(i);
    helper(q);
    subplot(1, 3, i); show(q);
end

function helper(q)
    confDump = sprintf('./output/%04.4f-%.4f-%d-%d.dump', ...
        q.beta, q.h, q.N, q.dimension);
    corrDump = 'noDump';
    saveFileName = '/dev/null';
    command = sprintf('../main --outputInt %d --skip %d --maxItr %d --J %.2f --h %.2f --beta %.3f --N %d --dimension %d --corrDump %s --confDump %s 2>%s', ...
        q.outputInt, q.preWait, (q.sample+q.preWait)*q.outputInt, ...
        q.J, q.h, q.beta, q.N, q.dimension, ...
        corrDump, confDump, saveFileName);
    
    tic
    system(command);
    toc
end

function show(q)
    confDump = sprintf('./output/%04.4f-%.4f-%d-%d.dump', ...
        q.beta, q.h, q.N, q.dimension);
    mat = load(confDump, '-ascii');
    imagesc(mat);
    title(['\beta=' num2str(q.beta)])
end

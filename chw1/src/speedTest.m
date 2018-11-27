clear; %clc

q.sample = 9000;
q.preWait = 1000;
q.outputInt = 1e4;
q.N = 14;
q.dimension = 3;

q.J = 1.0;
q.h = 0.0;

q.beta = 0.226;
helper(q);

function helper(q)
    confDump = 'noDump';
    corrDump = 'noDump';
    saveFileName = '/dev/null';
    command = sprintf('./main --outputInt %d --skip %d --maxItr %d --J %.2f --h %.2f --beta %.3f --N %d --dimension %d --corrDump %s --confDump %s 2>%s', ...
        q.outputInt, q.preWait, (q.sample+q.preWait)*q.outputInt, ...
        q.J, q.h, q.beta, q.N, q.dimension, ...
        corrDump, confDump, saveFileName);
    
    tic
    system(command);
    toc
end
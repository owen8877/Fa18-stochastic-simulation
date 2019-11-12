function timesMean = runTestbcc(qn, query)
    addpath lib
    
    system('make main_bcc');
    N = query.N;
    dimension = query.dimension;
    J = query.J;
    h = query.h;
    sample = query.sample;
    preWait = query.preWait;
    outputInt = query.outputInt;
    allTries = query.allTries;
    beta = query.beta;
    noDump = default(query, 'noDump', false);
    
    times = zeros(allTries, 1);
    time = 0;
    for j = 1:allTries
        clc;
        fprintf('** %.3f x%d (last:%.2fs) **\n', beta, j, time)

        if noDump
            corrDump = 'noDump';
        else
            corrDump = sprintf('./output/%s/%04.4f-%.4f-%d-%d.dump.%d', ...
                qn, beta, h, N, dimension, j);
        end
        saveFileName = sprintf('./output/%s/%04.4f-%.4f-%d-%d.mat.%d', ...
            qn, beta, h, N, dimension, j);
        command = sprintf('./main_bcc --outputInt %d --skip %d --maxItr %d --J %.2f --h %.2f --beta %.3f --N %d --dimension %d --corrDump %s 2>%s', ...
            outputInt, preWait, (sample+preWait)*outputInt, ...
            J, h, beta, N, dimension, ...
            corrDump, saveFileName);
        tic
        system(command);
        time = toc;
        times(j) = time;
    end
    timesMean = mean(times);
end


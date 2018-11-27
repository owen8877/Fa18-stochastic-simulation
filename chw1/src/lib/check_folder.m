function check_folder(path)
    if ~exist(path, 'dir')
        mkdir(path)
    end
end


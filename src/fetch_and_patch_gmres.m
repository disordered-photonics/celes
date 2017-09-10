function fetch_and_patch_gmres()

    out_dir = fullfile(tempdir, 'celes');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    in_dir = fullfile(matlabroot, 'toolbox', 'distcomp', 'array', '+parallel', '+internal', '+flowthrough');
    copyfile(fullfile(in_dir, 'private', 'iterapp.m'), out_dir);
    fileattrib(fullfile(out_dir, 'iterapp.m'), '+w');  % make it writable so that it can be overwritten or deleted

    copyfile(fullfile(in_dir, 'itermsg.m'), out_dir);
    fileattrib(fullfile(out_dir, 'itermsg.m'), '+w');

    path = fullfile(in_dir, 'gmres.m');
    v = fileread(path);  % read file into vector of chars
    newline_idxs = find(v == newline);  % find indices of newlines

    first_part_end_idx = newline_idxs(89);
    second_part_start_idx = newline_idxs(393);

    v_new = [...
        v(1:44), '_monitor', ...
        v(45:first_part_end_idx), newline, ...
        'msg = []; iternum = 0;', newline, ...
        v(first_part_end_idx:second_part_start_idx), ...
        '        iternum = iternum + 1; msg = convergence_message(msg, iternum, normr/n2b);', newline, ...
        v(second_part_start_idx:end)
    ];

    path = fullfile(out_dir, 'gmres_monitor.m');
    fd = fopen(path, 'w');
    fwrite(fd, v_new);
    fclose(fd);

    addpath(genpath(out_dir));

end

%  Copyright (c) 2017, Amos Egel (KIT), Lorenzo Pattelli (LENS)
%                      Giacomo Mazzamuto (LENS)
%  All rights reserved.
%
%  Redistribution and use in source and binary forms, with or without
%  modification, are permitted provided that the following conditions are met:
%
%  * Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
%
%  * Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
%
%  * Neither the name of the copyright holder nor the names of its
%    contributors may be used to endorse or promote products derived from
%    this software without specific prior written permission.
%
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
%  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%  POSSIBILITY OF SUCH DAMAGE.

%===============================================================================
%> @brief Copies MATLAB's gmres.m function to a temporary folder and adds status
%>        output messages to allow for a live monitoring of the convergence
%===============================================================================
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

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
%> @brief Compile the mex cuda function coupling_matrix_multiply_CUDA.cu
%>
%> @param lmax (int): maximal degree used in the SVWF expansion
%>
%> @param verbose (boolean): be verbose in compilation output?
%===============================================================================
function cuda_compile(lmax, verbose)

if exist('lmax','var') == 1
    lmaxmsg = sprintf(' with lmax=%d.\n',lmax);
else
    lmax = 3; % default lmax
    lmaxmsg = sprintf(' with default lmax.\n');
end

SRC_FILE = 'coupling_matrix_multiply_CUDA.cu';
DEFINES = strcat('-DLMAX=',int2str(lmax));

fprintf(1,['compiling CUDA-code',lmaxmsg])

OUT_FILE = [SRC_FILE(1:end-3), '_lmax', int2str(lmax)];
script_dir = fileparts(mfilename('fullpath'));
SRC_FILE = fullfile(script_dir, SRC_FILE);

verb_flag = '';

if nargin > 1
    if verbose
        verb_flag = '-v';
    end
end

if verLessThan('matlab','8.4') % check if version is less than R2015b
    myArch = computer('arch'); pathToOpts = fullfile(matlabroot, ...
        'toolbox', 'distcomp', 'gpu', ...
        'extern', 'src', 'mex', myArch,'gcc',...
        ['mex_CUDA_' myArch '.xml']);
    copyfile(pathToOpts,script_dir,'f')
    
    mex(verb_flag, DEFINES, '-outdir', script_dir, SRC_FILE)
else
    mexcuda(verb_flag, DEFINES, '-outdir', script_dir, '-output', OUT_FILE, SRC_FILE)
end
end

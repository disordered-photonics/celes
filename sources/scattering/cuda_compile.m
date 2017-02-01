%======================================================================
%> @brief Compile the mex cuda function coupling_matrix_multiply_CUDA.cu
%>
%> @param       lmax (int): maximal degree used in the SVWF expansion
%>
%> @param       verbose (boolean): be verbose in compilation output? 
%======================================================================
function cuda_compile(lmax, verbose)
    SRC_FILE = 'coupling_matrix_multiply_CUDA.cu';
    script_dir = fileparts(mfilename('fullpath'));
    SRC_FILE = fullfile(script_dir, SRC_FILE);

    if exist('lmax','var') == 1
        fprintf(1,'Building with lmax=%d.\n',lmax)
        DEFINES = strcat('-DLMAX=',int2str(lmax));
    else
        fprintf(1,'Building with default lmax.\n')
        DEFINES = '';
    end

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
        mexcuda(verb_flag, DEFINES, '-outdir', script_dir, SRC_FILE)
    end
end

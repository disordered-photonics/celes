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
%> @brief Numerically evaluate the coupling matrix multiply with some
%>        coefficient vector. The calculation is done on the GPU
%>
%> This routine is essentially a wrapper for the cuda mex function
%> coupling_matrix_multiply_CUDA
%>
%> @param   simulation (celes_simulation)
%>
%> @param   x (float array): coefficient vector to be multiplied with
%>          the coupling matrix
%>
%> @retval  Wx (float gpuArray): W*x, where W is the SVWF translation operator
%===============================================================================
function [Wx] = coupling_matrix_multiply(simulation,x,varargin)

if isempty(varargin)
    verbose = false;
else
    verbose = varargin{1};
end

if verbose
    fprintf('prepare particle coupling ... ')
end
preparation_timer = tic;

lmax = simulation.numerics.lmax;
NS = simulation.input.particles.number;
PlmCoeffTable = gpuArray(single(simulation.numerics.Plm_coeff_table));

real_hTab = simulation.h3_table.real_h3;
imag_hTab = simulation.h3_table.imag_h3;
rResol = single(simulation.numerics.particleDistanceResolution);

real_x = gpuArray(single(real(x)));
imag_x = gpuArray(single(imag(x)));

real_ab5Tab = [];
imag_ab5Tab = [];

loopCounter = 1;
for tau1 = 1:2
    for l1 = 1:lmax
        for m1 = -l1:l1
            j1 = multi2single_index(1,tau1,l1,m1,lmax);
            for tau2 = 1:2
                for l2 = 1:lmax
                    for m2 = -l2:l2
                        j2 = multi2single_index(1,tau2,l2,m2,lmax);
                        for p = max(abs(m1-m2),abs(l1-l2)+abs(tau1-tau2)):(l1+l2)
                            real_ab5Tab(end+1) = single(real(simulation.tables.translationTable.ab5(j2,j1,p+1)));
                            imag_ab5Tab(end+1) = single(imag(simulation.tables.translationTable.ab5(j2,j1,p+1)));
                            loopCounter = loopCounter+1;
                        end
                    end
                end
            end
        end
    end
end

real_ab5Tab = gpuArray(single(real_ab5Tab));
imag_ab5Tab = gpuArray(single(imag_ab5Tab));

part_pos = gpuArray(single(transpose(simulation.input.particles.positionArray)));

preparation_time = toc(preparation_timer);
if verbose
    fprintf('done in %f seconds.\n', preparation_time)
    fprintf('compute particle coupling ... ')
end

coupling_timer = tic;

mexfilename = ['coupling_matrix_multiply_CUDA_lmax', int2str(lmax)];
coupling_matrix_multiply_CUDA = str2func(mexfilename);

[real_Wx,imag_Wx] = coupling_matrix_multiply_CUDA(real_x, imag_x, ...
                                                  real_hTab, imag_hTab, ...
                                                  PlmCoeffTable, ...
                                                  real_ab5Tab, imag_ab5Tab, ...
                                                  part_pos, int32(NS), rResol);

Wx = real_Wx + 1i*imag_Wx;

coupling_time = toc(coupling_timer);
if verbose
    fprintf('done in %f seconds.\n', coupling_time)
end
end

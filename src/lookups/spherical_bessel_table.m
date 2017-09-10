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

%======================================================================
%> @brief Generate tabulated data of the spherical Hankel function
%>
%> @param simulation (celes_simulation)
%>
%> @retval h3_table (gpuArray): h3_table.real(i,j) contains the real part
%> of the spherical Hankel function of first kind of order i-1, evaluated 
%> at simulation.input.k_medium * simulation.lookupParticleDistances,
%> h3_table.imag(i,j) accordingly.
%======================================================================
function h3_table = spherical_bessel_table(simulation)

h3_table.real_h3=gpuArray.zeros(2*simulation.numerics.lmax+1,length(simulation.lookupParticleDistances),'single');
h3_table.imag_h3=gpuArray.zeros(2*simulation.numerics.lmax+1,length(simulation.lookupParticleDistances),'single');
for p=0:2*simulation.numerics.lmax
    spbs = sph_bessel(3,p,simulation.input.k_medium * simulation.lookupParticleDistances);
    h3_table.real_h3(p+1,:) = real(spbs);
    h3_table.imag_h3(p+1,:) = imag(spbs);
end

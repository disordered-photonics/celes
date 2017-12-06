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
%> @brief Cartesian components of a vector plane wave with unit amplitude
%>
%> @param R (Nx3 float array): positions where to evaluate wave
%> @param k (float): wavenumber
%> @param alpha (float): azimuthal angle of propagation direction
%> @param beta (float): polar angle of propagation direction
%> @param pol (int): polarization (1=TE, 2=TM)
%>
%> @retval E_components (cell array): {Ex,Ey,Ez}
%===============================================================================
function [E_components] = PVWF_components(R,k,alpha,beta,pol)
% alpha, beta need to be row vectors
% R need to be Nx3

kvec = k*[sin(beta).*cos(alpha);sin(beta).*sin(alpha);cos(beta)];
E = exp(1i*R*kvec);
if pol==1  % TE
    % E_alpha = [-sin(alpha);cos(alpha);0];
    Ex = -sin(alpha).*E;
    Ey = cos(alpha).*E;
    Ez = zeros(size(E),'like',E);
else
    % E_beta = [cos(alpha)*kz/k;sin(alpha)*kz/k;-kpar/k];
    Ex = (cos(alpha).*cos(beta)).*E;
    Ey = (sin(alpha).*cos(beta)).*E;
    Ez = -sin(beta).*E;
end

E_components = {Ex,Ey,Ez};
end

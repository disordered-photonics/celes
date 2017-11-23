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
%> @brief Normalized associated Legendre function for complex arguments
%>
%> The algorithm to compute the functions was taken from
%> "Light Scattering by Systems of Particles, Null-Field Method with
%> Discrete Sources: Theory and Programs"
%> by A. Doicu, T. Wriedt, and Y.A. Eremin
%>
%> For the normalization convention, see the \ref theory section.
%>
%> @param ct (float array): cos(theta)
%> @param st (float array): sin(theta)
%> @param lmax (int): maximal degree l of P_l^m
%>
%> @retval plm (cell array): plm{l+1,m+1} contains P_l^m(cos(theta)).
%===============================================================================
function plm = legendre_normalized_trigon(ct,st,lmax)
%
% Return the associated Legendre functions of cos(theta)=kz/k
%
% input arguments:
% - ct: cos(theta)
% - st: sin(theta)
% - lmax:  maximal polar angle quantum number
%
% output:
% - plm:   cell array of dimension lmax+1 x lmax+1
%          plm{l+1,m+1} contains P_l^m(cos(theta))
%
% The algorithm to compute the functions is taken from
% "Light Scattering by Systems of Particles, Null-Field Method with
% Discrete Sources: Theory and Programs"
% by A. Doicu, T. Wriedt, and Y.A. Eremin

plm = cell(lmax+1,lmax+1);  % first index: l+1, second index: m+1, inside each cell: same dimension as kpar

plm{1,1} = ones(size(ct),'like',ct)*sqrt(2)/2;  % P_0
plm{2,1} = sqrt(3/2)*ct; % P_1

for l = 1:lmax-1  % l+1
    plm{l+2,1} = 1/(l+1)*sqrt((2*l+1)*(2*l+3))*ct.*plm{l+1,1}- ...
                 l/(l+1)*sqrt((2*l+3)/(2*l-1))*plm{l,1};
end

for m = 1:lmax
    plm{m,m+1} = zeros(size(ct),'like',ct);
    plm{m+1,m+1} = sqrt((2*m+1)/2/factorial(2*m))*prod((2*m-1):-2:1)*st.^m;
    for l = m:lmax-1
        plm{l+2,m+1} = sqrt((2*l+1)*(2*l+3)/(l+1-m)/(l+1+m))*ct.*plm{l+1,m+1}- ...
                         sqrt((2*l+3)*(l-m)*(l+m)/(2*l-1)/(l+1-m)/(l+1+m))*plm{l,m+1};
    end
end
end

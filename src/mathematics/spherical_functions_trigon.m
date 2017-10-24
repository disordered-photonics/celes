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
%> @brief Spherical functions pi and tau
%>
%> The algorithm to compute the functions was taken from
%> "Light Scattering by Systems of Particles, Null-Field Method with
%> Discrete Sources: Theory and Programs"
%> by A. Doicu, T. Wriedt, and Y.A. Eremin
%>
%> For the definition of the pi and tau functions, see the \ref theory
%> section
%>
%> @param ct (float array): cos(theta)
%> @param st (float array): sin(theta)
%> @param lmax (int): maximal degree l of P_l^m
%>
%> @retval pilm (cell array): pilm{l+1,m+1} contains pi_l^m(cos(theta))
%> @retval taulm (cell array): taulm{l+1,m+1} contains pi_l^m(cos(theta))
%===============================================================================
function [pilm,taulm] = spherical_functions_trigon(ct,st,lmax)

plm = cell(lmax+1,lmax+1);  % first index: l+1, second index: m+1, inside each cell: same dimension as kpar
pilm = cell(lmax+1,lmax+1);  % first index: l+1, second index: m+1, inside each cell: same dimension as kpar
taulm = cell(lmax+1,lmax+1);  % first index: l+1, second index: m+1, inside each cell: same dimension as kpar
pprimel0 = cell(lmax+1);  % first index: l+1, second index: m+1, inside each cell: same dimension as kpar

plm{1,1} = sqrt(2)/2*ones(size(ct),'like',ct);  % P_0
plm{2,1} = sqrt(3/2)*ct; % P_1

pilm{1,1} = zeros(size(ct),'like',ct);  % pi_0^0
pilm{2,1} = zeros(size(ct),'like',ct); % pi_1^0

pprimel0{1} = zeros(size(ct),'like',ct);  % P'_0^0
pprimel0{2} = sqrt(3)*plm{1,1}; % P'_1^0

taulm{1,1} = -st.*pprimel0{1};  % tau_0^0
taulm{2,1} = -st.*pprimel0{2};  % tau_1^0

for l = 1:lmax-1  % l+1
    lp1=l+1;  % as index for cell array
    plm{lp1+1,1} = 1/(l+1)*sqrt((2*l+1)*(2*l+3))*ct.*plm{lp1,1}- ...
                   l/lp1*sqrt((2*l+3)/(2*l-1))*plm{lp1-1,1};
    pilm{lp1+1,1} = zeros(size(ct),'like',ct);
    coeff = sqrt((2*(l+1)+1)/(2*(l+1)-1));
    pprimel0{lp1+1} = (l+1)*coeff*plm{lp1,1}+coeff*ct.*pprimel0{lp1};
    taulm{lp1+1,1} = -st.*pprimel0{lp1+1};
end

for m=1:lmax
    mp1=m+1;
    plm{mp1-1,mp1}=zeros(size(ct),'like',ct);
    pilm{mp1-1,mp1}=zeros(size(ct),'like',ct);
    coeff = sqrt((2*m+1)/2/factorial(2*m))*prod((2*m-1):-2:1);
    plm{mp1,mp1}=coeff*st.^m;
    pilm{mp1,mp1}=coeff*st.^(m-1);
    taulm{mp1,mp1}=m*ct.*pilm{mp1,mp1};
    for l=m:lmax-1
        lp1=l+1;
        coeff1 = sqrt((2*l+1)*(2*l+3)/(l+1-m)/(l+1+m))*ct;
        coeff2 = sqrt((2*l+3)*(l-m)*(l+m)/(2*l-1)/(l+1-m)/(l+1+m));
        plm{lp1+1,mp1}=coeff1.*plm{lp1,mp1} - coeff2*plm{lp1-1,mp1};
        pilm{lp1+1,mp1}=coeff1.*pilm{lp1,mp1} - coeff2*pilm{lp1-1,mp1};
        taulm{lp1+1,mp1}=(l+1)*ct.*pilm{lp1+1,mp1}- ...
           (l+1+m)*sqrt((2*(l+1)+1)*(l+1-m)/(2*(l+1)-1)/(l+1+m)).*pilm{lp1,mp1};
    end
end
end

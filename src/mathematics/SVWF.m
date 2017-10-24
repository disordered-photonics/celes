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
%> @brief Spherical vector wave functions
%>
%> The convention for the definition of the SVWFs is taken from
%> "Light Scattering by Systems of Particles, Null-Field Method with
%> Discrete Sources: Theory and Programs"
%> by A. Doicu, T. Wriedt, and Y.A. Eremin
%>
%> See also the \ref theory section.
%>
%> @param k (float): wavenumber
%> @param R (3xN float array): positions where to evaluate the field in the
%>        format [x1,...,xN; y1,...,yN; z1,...,zN]
%> @param nu (int): regular (nu=1) or outgoing (nu=3) SVWF
%> @param tau (int): polarization (spherical TE, i.e., rxM=0 is tau=1,
%>        spherical TM is tau=2)
%> @param l (int): polar quantum number (degree), l=1,...
%> @param m (int): azimuthal quantum number (order), m=-l,...,l
%>
%> @retval N (float array): field components in format
%>         [Ex1,...,Exn; Ey1,...,Eyn; Ez1,...,Ezn]
%===============================================================================
function N = SVWF(k,R,nu,tau,l,m)
% R in the format [x1,...,xn; y1,...,yn; z1,...,zn]

% spherical coordinates
r = sqrt(sum(R.^2,2));
e_r = R./r;

ct = e_r(:,3);
st = sqrt(1-ct.^2);
phi = atan2(R(:,2),R(:,1));

e_theta = [ct.*cos(phi),ct.*sin(phi),-st];
e_phi = [-sin(phi),cos(phi),zeros(size(phi),'like',phi)];

% spherical functions
[p_all] = legendre_normalized_trigon(ct,st,l);
[pi_all,tau_all] = spherical_functions_trigon(ct,st,l);

P_lm = p_all{l+1,abs(m)+1};

pi_lm = pi_all{l+1,abs(m)+1};

tau_lm = tau_all{l+1,abs(m)+1};

z=sph_bessel(nu,l,k*r);

dxxz=dx_xz(nu,l,k*r);

eimphi=exp(1i*m*phi);

kr = k*r;

% SVWFs
if tau == 1 %select M
    N = 1/sqrt(2*l*(l+1)) * z .* (1i*m*pi_lm.*e_theta-tau_lm.*e_phi) .* eimphi;
else %select N
    N = 1/sqrt(2*l*(l+1)) * (l*(l+1)*z./kr.*P_lm.*e_r+ ...
        dxxz./kr.* (tau_lm.*e_theta + 1i*m*pi_lm.*e_phi)) .* eimphi;
end
end

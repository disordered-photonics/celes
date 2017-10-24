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
%> @brief Evaluate the Mie coefficients, see Bohren and Huffman:
%>        "Absorption and Scattering of Light by Small Particles", equations
%>        (4.52) and (4.53).
%>
%> @param   tau (int): SVWF polarization (1=TE, 2=TM)
%>
%> @param   l (int): SVWF degree (polar quantum number)
%>
%> @param   km (float): wavenumber in surrounding medium
%>
%> @param   kS (float): wavenumber inside particle
%>
%> @param   R (float): radius of sphere
%>
%> @retval  Q (float): Mie coefficient
%===============================================================================
function Q = T_entry(tau,l,kM,kS,R,varargin)

% The conventions are according to Bohren and Huffman's textbook.
% There is a twist: for tau=1 we have Q=-b and for tau=2 we have Q=-a in
% the case of scattered field Mie coefficients.
% This has to do with the definition of a and b in BH, see also
% Mishchenko: "Scattering, Absorption and Emission of Light by Small
% Paritcles", equations (5.42) and (5.43).

m = kS/kM;
x = kM*R;
mx = kS*R;

jx = sph_bessel(1,l,x);
jmx = sph_bessel(1,l,mx);
hx = sph_bessel(3,l,x);
djx = dx_xz(1,l,x);
djmx = dx_xz(1,l,mx);
dhx = dx_xz(3,l,x);

if(isempty(varargin))
    varargin={'scattered'};
end

switch lower(varargin{1})
    case 'scattered'
        if tau==1
            Q = -(jmx*djx-jx*djmx)/(jmx*dhx-hx*djmx); % -b
        else
            Q = -(m^2*jmx*djx-jx*djmx)/(m^2*jmx*dhx-hx*djmx); % -a
        end
    case 'internal'
        if tau==1
            Q = (jx*dhx-hx*djx)/(jmx*dhx-hx*djmx); % c
        else
            Q = (m*jx*dhx-m*hx*djx)/(m^2*jmx*dhx-hx*djmx); % d
        end
end
end

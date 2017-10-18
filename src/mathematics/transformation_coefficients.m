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
%> @brief Transformation operator B to map from spherical to plane vector wave
%>        functions and vice versa
%>
%> For the definition of B, see the \ref theory section.
%>
%> @param       pilm (cell array): angular function pi as returned by
%>                  spherical_functions_trigon()
%>
%> @param       taulm (cell array): angular function tau as returned by
%>                  spherical_functions_trigon()
%>
%> @param       tau (int): SVWF polarization (1=TE, 2=TM)
%>
%> @param       l (int): polar quantum number (degree), l=1,...
%>
%> @param       m (int): azimuthal quantum number (order), m=-l,...,l
%>
%> @param       pol (int): PVWF polarization (1=TE, 2=TM)
%>
%> @param       dagkey (string): Optional: keyword 'dagger' to compute B^dagger
%>
%> @retval      B (float array): B operator, same dimension as entries of pilm
%======================================================================
function B = transformation_coefficients(pilm,taulm,tau,l,m,pol,varargin)
% Transformation matrix B and B^\dagger to transform plane into spherical vector
% wave functions and vice versa.

if isempty(varargin)
    ifac = 1i;
elseif strcmpi(varargin{1},'dagger')
    ifac = -1i;
else
    error('dagger or not?')
end

if tau == pol
    spher_fun = taulm{l+1,abs(m)+1};
else
    spher_fun = m*pilm{l+1,abs(m)+1};
end

B = -1/(ifac)^(l+1)/sqrt(2*l*(l+1))*(ifac*(pol==1)+(pol==2))*spher_fun;
end

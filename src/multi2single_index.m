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
%> @brief Map the SVWF coefficients and the sphere number to the multi-index
%>        jmult
%>
%> @param       jS (int): Particle number
%> @param       tau (int): SVWF polarization
%> @param       l (int): SVWF degree (polar quantum number)
%> @param       m (int): SVWF order (azimuthal quantum number)
%> @param       lmax (int): SVWF degree cutoff
%> @retval      jmult(int): multi-index
%===============================================================================

function jmult = multi2single_index(jS,tau,l,m,lmax)
% jmult = multi2single_index(jS,tau,l,m,lmax)
%
% Return one index subsuming over the indices:
% - jS=1,2,...  (sphere number)
% - tau=1,2     (polarization: 1=TE commonly denoted M, 2=TM commonly denoted N)
% - l=1,..,lmax (polar eigenvalue)
% - m=-l,..,l   (azimuthal eigenvalue)
%
% The multiindex is counted according to the following scheme:
% Count first m, then l, then tau, then jS.
%
%     jmult   |   jS,tau,l,m
%    -----------------------
%         1   |   1,1,1,-1
%         2   |   1,1,1,0
%         3   |   1,1,1,1
%         4   |   1,1,2,-2
%         5   |   1,1,2,-1
%         .   |    ....
%         .   |   1,1,lmax,lmax
%         .   |   1,2,1,-1
%         .   |   1,2,1,0
%         .   |    ....
%         .   |   1,2,lmax,lmax
%         .   |   2,1,1,-1
%         .   |    ....
%   jmult_max |   1,2,1,0

jmult = (jS-1)*2*lmax*(lmax+2)+(tau-1)*lmax*(lmax+2)+(l-1)*(l+1)+m+l+1;
end

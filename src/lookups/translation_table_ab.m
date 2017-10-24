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
%> @brief Generate tabulated coefficients for the SVWF translation operator
%>
%> @param lmax (int): truncation degree of the expansion in SVWF
%>
%> @retval translation (struct): structure containing the a5 and b5 coefficients
%>         for the SVWF addition theorem, see \ref theory.
%===============================================================================
function translation = translation_table_ab(lmax)

jmax = jmult_max(1,lmax);

translation.ab5 = zeros(jmax,jmax,2*lmax+1,'single');

for tau1 = 1:2
    for l1 = 1:lmax
        for m1 = -l1:l1
            j1 = multi2single_index(1,tau1,l1,m1,lmax);
            for tau2 = 1:2
                for l2 = 1:lmax
                    for m2 = -l2:l2
                        j2 = multi2single_index(1,tau2,l2,m2,lmax);
                        for p = 0:2*lmax
                            if tau1 == tau2
                                translation.ab5(j1,j2,p+1) = (1i)^(abs(m1-m2)-abs(m1)-abs(m2)+l2-l1+p) * (-1)^(m1-m2) ...
                                    * sqrt((2*l1+1)*(2*l2+1)/(2*l1*(l1+1)*l2*(l2+1))) * (l1*(l1+1)+l2*(l2+1)-p*(p+1)) * sqrt(2*p+1) ...
                                    * Wigner3j([l1,l2,p],[m1,-m2,-m1+m2]) * Wigner3j([l1,l2,p],[0,0,0]);
                            elseif p > 0
                                translation.ab5(j1,j2,p+1) = (1i)^(abs(m1-m2)-abs(m1)-abs(m2)+l2-l1+p) * (-1)^(m1-m2) ...
                                    * sqrt((2*l1+1)*(2*l2+1)/(2*l1*(l1+1)*l2*(l2+1))) * sqrt((l1+l2+1+p)*(l1+l2+1-p)*(p+l1-l2)*(p-l1+l2)*(2*p+1)) ...
                                    * Wigner3j([l1,l2,p],[m1,-m2,-m1+m2]) * Wigner3j([l1,l2,p-1],[0,0,0]);
                            end
                        end
                    end
                end
            end
        end
    end
end
end

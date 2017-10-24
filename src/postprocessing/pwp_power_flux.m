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
%> @brief Compute the power flux for a field specified by a certain plane
%>        wave pattern
%>
%> The power flux is evaluated according to the formula
%> \f$ P = \frac{2\pi^2}{\omega k \mu_0}\sum_{j=1}^2
%> \int\dd{\alpha}\int\dd{\beta}\sin(\beta) \abs{g_j(\alpha,\beta)}^2 \f$,
%> see the \ref theory section.
%>
%> @param   pwp (celes_planeWavePattern): plane wave pattern of the
%>          field for which the power flux shall be evaluated
%>
%> @param   simulation (celes_simulation object): simulation object
%>          that stores all parameters of the simulation
%>
%> @param   Optional: direction (string): select 'forward' or
%>          'backward'. If specified, only partial waves propagating
%>          in the forward or backward z-direction are considered
%>
%> @retval  P (float): power flux
%===============================================================================
function P=pwp_power_flux(pwp,simulation,varargin)


if isempty(varargin)
    beta = pwp.polarAngles;
    alpha = pwp.azimuthalAngles;
    g = pwp.expansionCoefficients;
else
    switch lower(varargin{1})
        case 'forward'
            forward_idcs = (cos(pwp.polarAngles) >= 0);
            beta = pwp.polarAngles(forward_idcs);
            alpha = pwp.azimuthalAngles;
            g = pwp.expansionCoefficients(:, forward_idcs);
        case 'backward'
            backward_idcs = (cos(pwp.polarAngles) <= 0);
            beta = pwp.polarAngles(backward_idcs);
            alpha = pwp.azimuthalAngles;
            g = pwp.expansionCoefficients(:, backward_idcs);
    end
end

bintgrnd = trapz(alpha,abs(g).^2);
intgrl = trapz(beta, bintgrnd.*sin(beta));
P = 2*pi^2/(simulation.input.omega*simulation.input.k_medium)*real(intgrl);
end

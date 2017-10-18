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
%> @brief Evaluate the electrical field of a Gaussian beam in the paraxial
%>        approximation. The expressions are taken from Wikipedia.
%>
%> @param     simulation (celes_simulation object): the Gaussian beam
%>            parameters are taken from simulation.input.initialField
%>
%> @retval    E (Nx3 float array): electric field in the format [Ex;Ey;Ez]
%>            Each column correspond to one field point specified in
%>            simulation.output.fieldPoints
%===============================================================================
function E = gaussian_beam_paraxial_approximation(simulation)

% relative coordinates
xG = simulation.input.initialField.focalPoint(1);
yG = simulation.input.initialField.focalPoint(2);
zG = simulation.input.initialField.focalPoint(3);
z = simulation.output.fieldPoints(:,3) - zG;
z(z==0)=1e-30;
rho = sqrt((simulation.output.fieldPoints(:,1)-xG).^2 + ...
           (simulation.output.fieldPoints(:,2)-yG).^2);

% wavelength
k = simulation.input.k_medium;
wl = 2*pi/k;

switch lower(simulation.input.initialField.polarization)
    case 'te'
        alphaG = simulation.input.initialField.azimuthalAngle;
    case 'tm'
        alphaG = simulation.input.initialField.azimuthalAngle-pi/2;
end
E0 = simulation.input.initialField.amplitude * ... % complex beam amplitude
    [-sin(alphaG), cos(alphaG),0];
w0 = simulation.input.initialField.beamWidth; % beam waist
zR = pi*w0^2/wl; % Rayleigh range
wz = w0 * sqrt(1+(z/zR).^2); % spot size parameter
Rz = z.*(1+(zR./z).^2); % radius of curvature
psiz = atan(z/zR); % Gouy phase

E = sign(cos(simulation.input.initialField.polarAngle))* ...
    E0 * w0./wz .* exp(-rho.^2./wz.^2).* ...
    exp(1i*(k*z + k*rho.^2./(2*Rz) - psiz));
end

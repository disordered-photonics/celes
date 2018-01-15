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
%> @brief Expansion coefficients of a plane wave under normal incidence as the
%>        initial field in terms of regular spherical vector wave functions
%>        relative to the particles centers
%>
%> @param simulation (celes_simulation)
%>
%> @retval aI (device_array (CPU or GPU) of dimension NS x nmax, single
%>         precision)
%===============================================================================
function aI = initial_field_coefficients_planewave(simulation)

lmax = simulation.numerics.lmax;
E0 = simulation.input.initialField.amplitude;
k = simulation.input.k_medium;

beta = simulation.input.initialField.polarAngle;
cb = cos(beta);
sb = sin(beta);
alpha = simulation.input.initialField.azimuthalAngle;

% pi and tau symbols for transformation matrix B_dagger
[pilm,taulm] = spherical_functions_trigon(cb,sb,lmax);  % Nk x 1

% cylindrical coordinates for relative particle positions
relativeParticlePositions = simulation.input.particles.positionArray - simulation.input.initialField.focalPoint;
kvec = k*[sb*cos(alpha);sb*sin(alpha);cb];
eikr = exp(1i*relativeParticlePositions*kvec);

clear k beta cb sb kvec relativeParticlePositions % clean up some memory?

% compute initial field coefficients
aI = simulation.numerics.deviceArray(zeros(simulation.input.particles.number,simulation.numerics.nmax,'single'));
for m=-lmax:lmax
    for tau=1:2
        for l=max(1,abs(m)):lmax
            n=multi2single_index(1,tau,l,m,lmax);
            aI(:,n) = 4 * E0 * exp(-1i*m*alpha) .* eikr .* transformation_coefficients(pilm,taulm,tau,l,m,simulation.input.initialField.pol,'dagger');
        end
    end
end
end

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
%> @brief Expansion coefficients of a Gaussian beam under normal incidence as
%>        the initial field in terms of regular spherical vector wave functions
%>        relative to the particles centers
%>
%> @param simulation (celes_simulation)
%>
%> @retval aI (device_array (CPU or GPU) of dimension NS x nmax, single
%>         precision)
%===============================================================================
function aI = initial_field_coefficients_wavebundle_normal_incidence(simulation)

%--------------------------------------------------------------------------
% Cheatsheet for the dimensions of the arrays used in this function:
% a ~ NS x nmax
% a(:,n) = (eikz.*J) * (B.*gauss)
% B ~ Nk x 1
% k ~ Nk x 1
% eimph ~ NS x 1
% eikz ~ NS x Nk
% J ~ NS x Nk
% order of indices: jS, jk, jn
%--------------------------------------------------------------------------

lmax=simulation.numerics.lmax;
E0 = simulation.input.initialField.amplitude;
k = simulation.input.k_medium;
w = simulation.input.initialField.beamWidth;
prefac = E0*k^2*w^2/pi;
switch lower(simulation.input.initialField.polarization)
    case 'te'
        alphaG = simulation.input.initialField.azimuthalAngle;
    case 'tm'
        alphaG = simulation.input.initialField.azimuthalAngle-pi/2;
end

% grid for integral
fullBetaArray = simulation.numerics.polarAnglesArray(:);
directionIdcs = ( sign(cos(fullBetaArray)) == sign(cos(simulation.input.initialField.polarAngle)) );
betaArray = fullBetaArray(directionIdcs);
dBeta = mean(diff(betaArray));
cb = cos(betaArray);
sb = sin(betaArray);

% gauss factor
gaussfac = exp(-w^2/4*k^2*sb.^2);   % Nk x 1
gaussfacSincos = gaussfac.*cb.*sb;

% pi and tau symbols for transformation matrix B_dagger
[pilm,taulm] = spherical_functions_trigon(cb,sb,lmax);  % Nk x 1

% cylindrical coordinates for relative particle positions
relativeParticlePositions = simulation.input.particles.positionArray - simulation.input.initialField.focalPoint;
rhoGi = sqrt(relativeParticlePositions(:,1).^2+relativeParticlePositions(:,2).^2); % NS x 1
phiGi = atan2(relativeParticlePositions(:,2),relativeParticlePositions(:,1)); % NS x 1
zGi = relativeParticlePositions(:,3); % NS x 1

clear fullBetaArray directionIdcs betaArray gaussfac relativeParticlePositions % clean up some memory?

% compute initial field coefficients
aI = simulation.numerics.deviceArray(zeros(simulation.input.particles.number,simulation.numerics.nmax,'single'));
for m=-lmax:lmax

    % calculate terms on the fly to avoid storing temporary NSxNk-sized variables
    % inefficient performance-wise, for some calculations are now repeated twice
    eikzI1 = pi*( ...
              exp(-1i*alphaG) * 1i^abs(m-1) * ( exp(-1i*(m-1)*phiGi) .* (exp(1i*zGi*k*cb.') .* besselj(abs(m-1),rhoGi*k*sb.')) ) ...
            + exp( 1i*alphaG) * 1i^abs(m+1) * ( exp(-1i*(m+1)*phiGi) .* (exp(1i*zGi*k*cb.') .* besselj(abs(m+1),rhoGi*k*sb.')) ) ...
                );

    eikzI2 = pi*1i*( ...
            - exp(-1i*alphaG) * 1i^abs(m-1) * ( exp(-1i*(m-1)*phiGi) .* (exp(1i*zGi*k*cb.') .* besselj(abs(m-1),rhoGi*k*sb.')) ) ...
            + exp( 1i*alphaG) * 1i^abs(m+1) * ( exp(-1i*(m+1)*phiGi) .* (exp(1i*zGi*k*cb.') .* besselj(abs(m+1),rhoGi*k*sb.')) ) ...
                    );
    for tau=1:2
        for l=max(1,abs(m)):lmax
            n=multi2single_index(1,tau,l,m,lmax);

            gaussSincosBDag1 = gaussfacSincos .* transformation_coefficients(pilm,taulm,tau,l,m,1,'dagger'); % Nk x 1
            gaussSincosBDag2 = gaussfacSincos .* transformation_coefficients(pilm,taulm,tau,l,m,2,'dagger'); % Nk x 1

            aI(:,n) = prefac * (...
                                eikzI1(:,2:end) * gaussSincosBDag1(2:end) + ...
                                eikzI1(:,1:(end-1)) * gaussSincosBDag1(1:end-1) + ...
                                eikzI2(:,2:end) * gaussSincosBDag2(2:end) + ...
                                eikzI2(:,1:(end-1)) * gaussSincosBDag2(1:end-1) ...
                                ) * dBeta / 2;
        end
    end
end
end

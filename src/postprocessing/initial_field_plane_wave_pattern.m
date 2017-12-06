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
%> @brief Compute the expansion of the initial field in plane vector wave
%>        functions
%>
%> @param   simulation (celes_simulation object): simulation for
%>          which the initial field pwp shall be computed
%>
%> @retval  pwp (1x2 cell array): pwp{1} and pwp{2} each contain a
%>          celes_planeWavePattern for the initial field in TE and in
%>          TM polarization, respectively
%===============================================================================
function pwp = initial_field_plane_wave_pattern(simulation)

switch simulation.input.initialField.type
    case 'gaussian beam'
        if simulation.input.initialField.normalIncidence

            initF = simulation.input.initialField;

            betaArray = simulation.numerics.polarAnglesArray;
            alphaArray = simulation.numerics.azimuthalAnglesArray;

            RG = simulation.input.initialField.focalPoint;
            E0 = initF.amplitude;
            w = initF.beamWidth;
            k = simulation.input.k_medium;

            for pol = 2:-1:1
            pwp{pol} = celes_planeWavePattern('polarAngles', betaArray, ...
                                              'azimuthalAngles', alphaArray, ...
                                              'k', simulation.input.k_medium, ...
                                              'expansionCoefficients', simulation.numerics.deviceArray(zeros(length(alphaArray),length(betaArray),'single')) ... % dima x dimb
                                              );
            end

            kxGrid = pwp{1}.kxGrid;
            kyGrid = pwp{1}.kyGrid;
            kzGrid = pwp{1}.kzGrid;

            emnikrg = exp(-1i*(kxGrid*RG(1) + kyGrid*RG(2) + kzGrid*RG(3)));  % Nalpha x Nbeta
            prefacCosBetGaussfac = E0*k^2*w^2/(4*pi)* cos(betaArray(:).') .* exp(-w^2/4*k^2*sin(betaArray(:).').^2)  .*  ( sign(cos(betaArray(:).')) == sign(cos(initF.polarAngle)) ) ;  % 1 x Nbeta, make sure only the correct direction contributes

            eikrgPrefacCosBetGaussfac = emnikrg.*prefacCosBetGaussfac; % Nalpha x Nbeta

            switch lower(simulation.input.initialField.polarization)
                case 'te'
                    alphaG = simulation.input.initialField.azimuthalAngle;
                case 'tm'
                    alphaG = simulation.input.initialField.azimuthalAngle-pi/2;
            end
            
            pwp{1}.expansionCoefficients = cos(alphaArray(:)-alphaG).*eikrgPrefacCosBetGaussfac;
            pwp{2}.expansionCoefficients = sign(cos(initF.polarAngle)) * (sin(alphaArray(:)-alphaG).*eikrgPrefacCosBetGaussfac);

        else
            error('this is not implemented')
        end
    otherwise
        error('this is not implemented')
end
end

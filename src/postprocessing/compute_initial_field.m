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
%> @brief Evaluate the near field of the initial excitation
%>
%> @param   simulation (celes_simulation object): simulation for
%>          which the initial field is to be computed
%>
%> @retval  E, H (Nx3 float arrays): electric and magnetic fields in the
%>          format [Ex;Ey;Ez]. Each column correspond to one field point
%>          specified in simulation.output.fieldPoints
%===============================================================================
function [E, H] = compute_initial_field(simulation)

E = zeros(size(simulation.output.fieldPoints),'single');  % Nx3
H = zeros(size(simulation.output.fieldPoints),'single');  % Nx3
numberOfFieldPoints = length(E(:,1));
R = simulation.numerics.deviceArray(single(simulation.output.fieldPoints));
k = simulation.input.k_medium;
n = simulation.input.mediumRefractiveIndex;

switch simulation.input.initialField.type
    case 'gaussian beam'

        if isfinite(simulation.input.initialField.beamWidth) && ...
                    simulation.input.initialField.beamWidth

            if simulation.output.initialFieldUseParaxialApproximation
                [E, H] = gaussian_beam_paraxial_approximation(simulation);
            else
                pwp = initial_field_plane_wave_pattern(simulation);
                for pol = 1:2

                    axIntegrandE = simulation.numerics.deviceArray(zeros(numberOfFieldPoints,length(pwp{pol}.azimuthalAngles),'single'));
                    ayIntegrandE = simulation.numerics.deviceArray(zeros(numberOfFieldPoints,length(pwp{pol}.azimuthalAngles),'single'));
                    azIntegrandE = simulation.numerics.deviceArray(zeros(numberOfFieldPoints,length(pwp{pol}.azimuthalAngles),'single'));
                    axIntegrandH = simulation.numerics.deviceArray(zeros(numberOfFieldPoints,length(pwp{pol}.azimuthalAngles),'single'));
                    ayIntegrandH = simulation.numerics.deviceArray(zeros(numberOfFieldPoints,length(pwp{pol}.azimuthalAngles),'single'));
                    azIntegrandH = simulation.numerics.deviceArray(zeros(numberOfFieldPoints,length(pwp{pol}.azimuthalAngles),'single'));

                    betArr = simulation.numerics.deviceArray(transpose(single(pwp{pol}.polarAngles(:))));
                    alphArr = simulation.numerics.deviceArray(transpose(single(pwp{pol}.azimuthalAngles(:))));
                    sinBet = sin(betArr);

                    for ja=1:length(alphArr)
                        %ja
                        [E_comps, H_comps] = PVWF_components(R,k,n,alphArr(ja),betArr,pol);

                        coeffs = sinBet.*pwp{pol}.expansionCoefficients(ja,:);
                        betaxIntegrand = coeffs.*E_comps{1};
                        betayIntegrand = coeffs.*E_comps{2};
                        betazIntegrand = coeffs.*E_comps{3};

                        axIntegrandE(:,ja) = trapz(betArr,betaxIntegrand,2);
                        ayIntegrandE(:,ja) = trapz(betArr,betayIntegrand,2);
                        azIntegrandE(:,ja) = trapz(betArr,betazIntegrand,2);

                        betaxIntegrand = coeffs.*H_comps{1};
                        betayIntegrand = coeffs.*H_comps{2};
                        betazIntegrand = coeffs.*H_comps{3};

                        axIntegrandH(:,ja) = trapz(betArr,betaxIntegrand,2);
                        ayIntegrandH(:,ja) = trapz(betArr,betayIntegrand,2);
                        azIntegrandH(:,ja) = trapz(betArr,betazIntegrand,2);
                    end
                    E(:,1) = E(:,1) + gather(trapz(alphArr,axIntegrandE,2));
                    E(:,2) = E(:,2) + gather(trapz(alphArr,ayIntegrandE,2));
                    E(:,3) = E(:,3) + gather(trapz(alphArr,azIntegrandE,2));
                    H(:,1) = H(:,1) + gather(trapz(alphArr,axIntegrandH,2));
                    H(:,2) = H(:,2) + gather(trapz(alphArr,ayIntegrandH,2));
                    H(:,3) = H(:,3) + gather(trapz(alphArr,azIntegrandH,2));
                end
            end
        end
    case 'plane wave'
        [E_comps, H_comps] = PVWF_components(R, k, n, ...
            simulation.input.initialField.azimuthalAngle, ...
            simulation.input.initialField.polarAngle, ...
            simulation.input.initialField.pol);
        
        E(:,1) = gather(E_comps{1});
        E(:,2) = gather(E_comps{2});
        E(:,3) = gather(E_comps{3});
        H(:,1) = gather(H_comps{1});
        H(:,2) = gather(H_comps{2});
        H(:,3) = gather(H_comps{3});
    otherwise
        error('this is not implemented')
end
end

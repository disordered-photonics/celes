%======================================================================
%> @brief Evaluate the near field of the initial excitation
%>
%> @param       simulation (celes_simulation object): simulation for
%>              which the initial field is to be computed
%>
%> @retval      E (Nx3 float array): electric field in the format [Ex;Ey;Ez]
%>              Each column correspond to one field point specified in 
%>              simulation.output.fieldPoints
%======================================================================
function E = compute_initial_field(simulation)

E = zeros(size(simulation.output.fieldPoints),'single');  % Nx3
numberOfFieldPoints = length(E(:,1));
R = simulation.numerics.deviceArray(single(simulation.output.fieldPoints));
k = simulation.input.k_medium;

switch simulation.input.initialField.type
    case 'gaussian beam'
        
        if isfinite(simulation.input.initialField.beamWidth) && simulation.input.initialField.beamWidth
            
            if simulation.output.initialFieldUseParaxialApproximation
                E = gaussian_beam_paraxial_approximation(simulation);
            else
                pwp = initial_field_plane_wave_pattern(simulation);
                for pol=1:2
                    
                    axIntegrand = simulation.numerics.deviceArray(zeros(numberOfFieldPoints,length(pwp{pol}.azimuthalAngles),'single'));
                    ayIntegrand = simulation.numerics.deviceArray(zeros(numberOfFieldPoints,length(pwp{pol}.azimuthalAngles),'single'));
                    azIntegrand = simulation.numerics.deviceArray(zeros(numberOfFieldPoints,length(pwp{pol}.azimuthalAngles),'single'));
                    
                    betArr = simulation.numerics.deviceArray(transpose(single(pwp{pol}.polarAngles(:))));
                    alphArr = simulation.numerics.deviceArray(transpose(single(pwp{pol}.azimuthalAngles(:))));
                    sinBet = sin(betArr);
                    
                    for ja=1:length(alphArr)
                        %ja
                        E_comps = PVWF_components(R,k,alphArr(ja),betArr,pol);
                        
                        coeffs = sinBet.*pwp{pol}.expansionCoefficients(ja,:);
                        betaxIntegrand = bsxfun(@times,coeffs,E_comps{1});
                        betayIntegrand = bsxfun(@times,coeffs,E_comps{2});
                        betazIntegrand = bsxfun(@times,coeffs,E_comps{3});
                        
                        axIntegrand(:,ja) = trapz(betArr,betaxIntegrand,2);
                        ayIntegrand(:,ja) = trapz(betArr,betayIntegrand,2);
                        azIntegrand(:,ja) = trapz(betArr,betazIntegrand,2);
                    end
                    E(:,1) = E(:,1) + gather(trapz(alphArr,axIntegrand,2));
                    E(:,2) = E(:,2) + gather(trapz(alphArr,ayIntegrand,2));
                    E(:,3) = E(:,3) + gather(trapz(alphArr,azIntegrand,2));
                end
            end
            
        else  % plane wave
            E_comps = PVWF_components(R,k,simulation.input.initialField.azimuthalAngle,simulation.input.initialField.polarAngle,simulation.input.initialField.pol);
            E(:,1) = gather(E_comps{1});
            E(:,2) = gather(E_comps{2});
            E(:,3) = gather(E_comps{3});
        end
        
    otherwise
        error('this is not implemented')
end
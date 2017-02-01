%======================================================================
%> @brief Compute the expansion of the initial field in plane vector wave
%> functions
%>
%> @param       simulation (celes_simulation object): simulation for
%>              which the initial field pwp shall be computed
%>
%> @retval      pwp (1x2 cell array): pwp{1} and pwp{2} each contain a
%>              celes_planeWavePattern for the initial field in TE and in
%>              TM polarization, respectively
%======================================================================
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
            
            pwp{1} = celes_planeWavePattern(betaArray,alphaArray,simulation.input.k_medium);
            pwp{2} = celes_planeWavePattern(betaArray,alphaArray,simulation.input.k_medium);
            
            kxGrid = pwp{1}.kxGrid;
            kyGrid = pwp{1}.kyGrid;
            kzGrid = pwp{1}.kzGrid;
            
            emnikrg = exp(-1i*(kxGrid*RG(1) + kyGrid*RG(2) + kzGrid*RG(3)));  % Nalpha x Nbeta
            prefacCosBetGaussfac = E0*k^2*w^2/(4*pi)* cos(betaArray(:).') .* exp(-w^2/4*k^2*sin(betaArray(:).').^2)  .*  ( sign(cos(betaArray(:).')) == sign(cos(initF.polarAngle)) ) ;  % 1 x Nbeta, make sure only the correct direction contributes
            
            eikrgPrefacCosBetGaussfac = bsxfun(@times,emnikrg,prefacCosBetGaussfac); % Nalpha x Nbeta
            
            switch simulation.input.initialField.polarization
                case 'TE'
                    alphaG = simulation.input.initialField.azimuthalAngle;
                case 'TM'
                    alphaG = simulation.input.initialField.azimuthalAngle-pi/2;
            end

            pwp{1}.expansionCoefficients = bsxfun(@times, cos(alphaArray(:)-alphaG), eikrgPrefacCosBetGaussfac);
            pwp{2}.expansionCoefficients = bsxfun(@times, sin(alphaArray(:)-alphaG), eikrgPrefacCosBetGaussfac);
            
        else
            error('this is not implemented')
        end
    otherwise
        error('this is not implemented')
end
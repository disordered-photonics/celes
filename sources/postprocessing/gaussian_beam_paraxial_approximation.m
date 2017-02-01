function E = gaussian_beam_paraxial_approximation(simulation)
%======================================================================
%> @brief Evaluate the electrical field of a Gaussian beam in the paraxial
%>        approximation. The expressions are taken from Wikipedia.
%>
%> @param     simulation (celes_simulation object): the Gaussian beam
%>            parameters are taken from simulation.input.initialField
%>
%> @retval    E (Nx3 float array): electric field in the format [Ex;Ey;Ez]
%>            Each column correspond to one field point specified in 
%>            simulation.output.fieldPoints
%======================================================================

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

switch simulation.input.initialField.polarization
    case 'TE'
        alphaG = simulation.input.initialField.azimuthalAngle;
    case 'TM'
        alphaG = simulation.input.initialField.azimuthalAngle-pi/2;
end
E0 = simulation.input.initialField.amplitude * ...  % complex beam amplitude
    [-sin(alphaG), cos(alphaG),0];
w0 = simulation.input.initialField.beamWidth; % beam waist
zR = pi*w0^2/wl; % Rayleigh range
wz = w0 * sqrt(1+(z/zR).^2); % spot size parameter
Rz = z.*(1+(zR./z).^2); % radius of curvature
psiz = atan(z/zR); % Gouy phase

E = E0 * w0./wz .* exp(-rho.^2./wz.^2) .* ...
    exp(1i*(k*z + k*rho.^2./(2*Rz) - psiz));
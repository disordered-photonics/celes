%======================================================================
%> @brief Power flux of incoming Gaussian beam under normal incidence
%>
%> @param simulation celes_simulation object
%>
%> @retval PI (float) incoming power
%======================================================================
function PI = initial_power_wavebundle_normal_incidence(simulation)

E0 = simulation.input.initialField.amplitude;
k = simulation.input.k_medium;
omega = simulation.input.omega;
w = simulation.input.initialField.beamWidth;
betaArray = simulation.numerics.polarAnglesArray;

integrand = sin(betaArray).*cos(betaArray).^2.*exp(-w^2/2*k^2*sin(betaArray).^2) .* (sign(cos(betaArray)) == sign(cos(simulation.input.initialField.polarAngle)));

PI = abs(E0)^2*pi*k^3*w^4/(4*omega) * trapz(betaArray,integrand);
%======================================================================
%> @brief Expansion coefficients of a plane wave under normal incidence 
%> as the initial field in terms of regular spherical vector wave functions
%> relative to the particles centers
%>
%> @param simulation (celes_simulation)
%>
%> @retval aI (device_array (cpu or gpu) of dimension NS x nmax, single
%> precision)
%======================================================================
function aI = initial_field_coefficients_planewave(simulation)

lmax=simulation.numerics.lmax;
E0 = simulation.input.initialField.amplitude;
k = simulation.input.k_medium;

beta = simulation.input.initialField.polarAngle;
alpha = simulation.input.initialField.azimuthalAngle;

% pi and tau symbols for transformation matrix B_dagger
[pilm,taulm] = spherical_functions_angular(beta,lmax);  % Nk x 1

% cylindrical coordinates for relative particle positions
relativeParticlePositions = bsxfun(@plus,simulation.input.particles.positionArray,-simulation.input.initialField.focalPoint);
kvec = k*[sin(beta)*cos(alpha);sin(beta)*sin(alpha);cos(beta)];
eikr = exp(1i*relativeParticlePositions*kvec);

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
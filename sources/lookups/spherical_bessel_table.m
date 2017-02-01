%======================================================================
%> @brief Generate tabulated data of the spherical Hankel function
%>
%> @param simulation (celes_simulation)
%>
%> @retval h3_table (gpuArray): h3_table.real(i,j) contains the real part
%> of the spherical Hankel function of first kind of order i-1, evaluated 
%> at simulation.input.k_medium * simulation.lookupParticleDistances,
%> h3_table.imag(i,j) accordingly.
%======================================================================
function h3_table = spherical_bessel_table(simulation)

h3_table.real_h3=gpuArray.zeros(2*simulation.numerics.lmax+1,length(simulation.lookupParticleDistances),'single');
h3_table.imag_h3=gpuArray.zeros(2*simulation.numerics.lmax+1,length(simulation.lookupParticleDistances),'single');
for p=0:2*simulation.numerics.lmax
    spbs = sph_bessel(3,p,simulation.input.k_medium * simulation.lookupParticleDistances);
    h3_table.real_h3(p+1,:) = real(spbs);
    h3_table.imag_h3(p+1,:) = imag(spbs);
end

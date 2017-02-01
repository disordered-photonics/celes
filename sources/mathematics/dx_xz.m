%======================================================================
%> @brief The derivative of x*j(x) where j is a spherical Bessel or Hankel
%> function
%>
%> @param nu (int): selects between spherical Bessel (nu=1) and Hankel (nu=3)
%> function
%> @param l (int): order of spherical Bessel/Hankel function
%> @param Z (complex float array): argument where the derivative be
%> evaluated
%>
%> @retval dxxz (complex float array): derivative of Z*sph_bessel(nu,l,Z)
%> with respect to Z
%======================================================================
function dxxz = dx_xz(nu,l,Z)
% dxxz = dx_xz(nu,l,Z)
% Return the derivative d_x (x*sph_bessel(nu,l,x))

% check 2014 02 06

dxxz = Z.*sph_bessel(nu,l-1,Z) - l*sph_bessel(nu,l,Z);
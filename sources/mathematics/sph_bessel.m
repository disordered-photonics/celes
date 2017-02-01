%======================================================================
%> @brief Spherical Bessel or Hankel function
%>
%> @param nu (int): selects between spherical Bessel (nu=1) and Hankel (nu=3)
%> function
%> @param l (int): order of spherical Bessel/Hankel function
%> @param Z (complex float array): argument 
%>
%> @retval sb (complex float array): spherical Bessel or Hankel function
%======================================================================
function sb = sph_bessel(nu,l,Z)
% function sb = sph_bessel(nu,l,Z)
% Spherical Bessel or Hankel function
%
% nu=1: spherical bessel function of first kind and order l
% nu=3: spherical hankel function of first kind and order l

% check 2014 02 06

if nu==1
    sb = sqrt(pi/2./Z).*besselj(l+1/2,Z); 
elseif nu==3
    sb = sqrt(pi/2./Z).*besselh(l+1/2,Z);
else
    error('nu must be 1 or 3')
end
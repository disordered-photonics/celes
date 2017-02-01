%======================================================================
%> @brief Spherical vector wave functions
%>
%> The convention for the definition of the SVWFs is taken from 
%> "Light Scattering by Systems of Particles, Null-Field Method with 
%> Discrete Sources: Theory and Programs" 
%> by A. Doicu, T. Wriedt, and Y.A. Eremin
%> 
%> See also the \ref theory section.
%> 
%> @param k (float): wavenumber
%> @param R (3xN float array): positions where to evaluate the field in the
%> format [x1,...,xN; y1,...,yN; z1,...,zN]
%> @param nu (int): regular (nu=1) or outgoing (nu=3) SVWF
%> @param tau (int): polarization (spherical TE, i.e., rxM=0 is tau=1,
%> spherical TM is tau=2)
%> @param l (int): polar quantum number (degree), l=1,...
%> @param m (int): azimuthal quantum number (order), m=-l,...,l
%>
%> @retval N (float array): field components in format 
%> [Ex1,...,Exn; Ey1,...,Eyn; Ez1,...,Ezn]
%======================================================================
function N = SVWF(k,R,nu,tau,l,m)
% R in the format [x1,...,xn; y1,...,yn; z1,...,zn]

% spherical coordinates
r = sqrt(R(:,1).^2+R(:,2).^2+R(:,3).^2);
theta = acos(R(:,3)./r);
phi = atan2(R(:,2),R(:,1));

e_r = R./[r,r,r];
e_theta = [cos(theta).*cos(phi),cos(theta).*sin(phi),-sin(theta)];
e_phi = [-sin(phi),cos(phi),r-r];

% spherical functions
[p_all] = legendre_normalized_angular(theta,l);
[pi_all,tau_all] = spherical_functions_angular(theta,l);

P_lm = p_all{l+1,abs(m)+1};
P_lm=[P_lm,P_lm,P_lm];

pi_lm = pi_all{l+1,abs(m)+1};
pi_lm=[pi_lm,pi_lm,pi_lm];

tau_lm = tau_all{l+1,abs(m)+1};
tau_lm=[tau_lm,tau_lm,tau_lm];

z=sph_bessel(nu,l,k*r);
z=[z,z,z];

dxxz=dx_xz(nu,l,k*r);
dxxz=[dxxz,dxxz,dxxz];

eimphi=exp(1i*m*phi);
eimphi=[eimphi,eimphi,eimphi];

kr = k*r;
kr = [kr,kr,kr];

% SVWFs
if tau==1  %select M
    N = 1/sqrt(2*l*(l+1)) * z .* (1i*m*pi_lm.*e_theta - tau_lm.*e_phi) .* eimphi;
else %select N
    N = 1/sqrt(2*l*(l+1)) * (l*(l+1)*z./kr.*P_lm.*e_r + dxxz./kr.* (tau_lm.*e_theta + 1i*m*pi_lm.*e_phi)) .* eimphi;
end


%======================================================================
%> @brief Cartesian components of a vector plane wave with unit amplitude
%>
%> @param R (Nx3 float array): positions where to evaluate wave
%> @param k (float): wavenumber
%> @param alpha (float): azimuthal angle of propagation direction
%> @param beta (float): polar angle of propagation direction
%> @param pol (int): polarization (1=TE, 2=TM)
%>
%> @retval E_components (cell array): {Ex,Ey,Ez}
%======================================================================
function [E_components] = PVWF_components(R,k,alpha,beta,pol)
% alpha, beta need to be row vectors
% R need to be Nx3

kvec = k*[sin(beta).*cos(alpha);sin(beta).*sin(alpha);cos(beta)];
E = exp(1i*R*kvec);
if pol==1  % TE
    % E_alpha = [-sin(alpha);cos(alpha);0];
    Ex = bsxfun(@times,-sin(alpha), E);
    Ey = bsxfun(@times,cos(alpha), E);
    Ez = E-E;
else
    % E_beta = [cos(alpha)*kz/k;sin(alpha)*kz/k;-kpar/k];
    Ex = bsxfun(@times,cos(alpha).*cos(beta), E);
    Ey = bsxfun(@times,sin(alpha).*cos(beta), E);
    Ez = bsxfun(@times,-sin(beta), E);
end

E_components = {Ex,Ey,Ez};
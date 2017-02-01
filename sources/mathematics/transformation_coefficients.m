%======================================================================
%> @brief Transformation operator B to map from spherical to plane vector
%> wave functions and vice versa
%>
%> For the definition of B, see the \ref theory section.
%>
%> @param       pilm (cell array): angular function pi as returned by
%>                  spherical_functions_angular()
%>
%> @param       taulm (cell array): angular function tau as returned by
%>                  spherical_functions_angular()
%>
%> @param       tau (int): SVWF polarization (1=TE, 2=TM)
%>
%> @param       l (int): polar quantum number (degree), l=1,...
%>
%> @param       m (int): azimuthal quantum number (order), m=-l,...,l
%>
%> @param       pol (int): PVWF polarization (1=TE, 2=TM)
%>
%> @param       dagkey (string): Optional: keyword 'dagger' to compute 
%>                  B^dagger
%>
%> @retval      B (float array): B operator, same dimension as entries of pilm
%======================================================================
function B = transformation_coefficients(pilm,taulm,tau,l,m,pol,varargin)
% Transformation matrix B and B^\dagger to transform plane into spherical vector wave
% functions and vice versa.

if isempty(varargin)
    ifac = 1i;
elseif strcmp(varargin{1},'dagger')
    ifac = -1i;
else
    error('dagger or not?')
end

if tau==pol
    spher_fun = taulm{l+1,abs(m)+1};
else
    spher_fun = m*pilm{l+1,abs(m)+1};
end

B = -1/(ifac)^(l+1)/sqrt(2*l*(l+1))*(ifac*(pol==1)+(pol==2))*spher_fun;
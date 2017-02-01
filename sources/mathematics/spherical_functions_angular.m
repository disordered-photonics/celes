%======================================================================
%> @brief Spherical functions pi and tau
%>
%> The algorithm to compute the functions was taken from 
%> "Light Scattering by Systems of Particles, Null-Field Method with 
%> Discrete Sources: Theory and Programs" 
%> by A. Doicu, T. Wriedt, and Y.A. Eremin
%> 
%> For the definition of the pi and tau functions, see the \ref theory
%> section
%> 
%> @param theta (float array): polar angles
%> @param lmax (int): maximal degree l of P_l^m
%>
%> @retval pilm (cell array): pilm{l+1,m+1} contains pi_l^m(cos(theta))
%> @retval taulm (cell array): taulm{l+1,m+1} contains pi_l^m(cos(theta))
%======================================================================
function [pilm,taulm] = spherical_functions_angular(theta,lmax)
% [pilm,taulm] = legendre_normalized(theta,lmax)
st = sin(theta);
ct = cos(theta);

plm = cell(lmax+1,lmax+1);  % first index: l+1, second index: m+1, inside each cell: same dimension as kpar
pilm = cell(lmax+1,lmax+1);  % first index: l+1, second index: m+1, inside each cell: same dimension as kpar
taulm = cell(lmax+1,lmax+1);  % first index: l+1, second index: m+1, inside each cell: same dimension as kpar
pprimel0 = cell(lmax+1);  % first index: l+1, second index: m+1, inside each cell: same dimension as kpar

plm{1,1}=theta-theta+sqrt(2)/2;  % P_0
plm{2,1} = sqrt(3/2)*ct; % P_1

pilm{1,1} = theta-theta;  % pi_0^0
pilm{2,1} = theta-theta; % pi_1^0

pprimel0{1} = theta-theta;  % P'_0^0
pprimel0{2} = sqrt(3)*plm{1,1}; % P'_1^0

taulm{1,1} = -st.*pprimel0{1};  % tau_0^0
taulm{2,1} = -st.*pprimel0{2};  % tau_1^0

for l=1:lmax-1  % l+1
    lp1=l+1;  % as index for cell array
    plm{lp1+1,1} = 1/(l+1)*sqrt((2*l+1)*(2*l+3))*ct.*plm{lp1,1} - l/lp1*sqrt((2*l+3)/(2*l-1))*plm{lp1-1,1};
    pilm{lp1+1,1} = theta-theta;
    pprimel0{lp1+1} = (l+1)*sqrt((2*(l+1)+1)/(2*(l+1)-1))*plm{lp1,1}+sqrt((2*(l+1)+1)/(2*(l+1)-1))*ct.*pprimel0{lp1};
    taulm{lp1+1,1} = -st.*pprimel0{lp1+1};    
end

for m=1:lmax
    mp1=m+1;
    plm{mp1-1,mp1}=theta-theta;
    pilm{mp1-1,mp1}=theta-theta;
    plm{mp1,mp1}=sqrt((2*m+1)/2/factorial(2*m))*factd2(2*m-1)*st.^m;
    pilm{mp1,mp1}=sqrt((2*m+1)/2/factorial(2*m))*factd2(2*m-1)*st.^(m-1);
    taulm{mp1,mp1}=m*ct.*pilm{mp1,mp1};
    for l=m:lmax-1
        lp1=l+1;
        plm{lp1+1,mp1}=sqrt((2*l+1)*(2*l+3)/(l+1-m)/(l+1+m))*ct.*plm{lp1,mp1} - sqrt((2*l+3)*(l-m)*(l+m)/(2*l-1)/(l+1-m)/(l+1+m))*plm{lp1-1,mp1};
        pilm{lp1+1,mp1}=sqrt((2*l+1)*(2*l+3)/(l+1-m)/(l+1+m))*ct.*pilm{lp1,mp1} - sqrt((2*l+3)*(l-m)*(l+m)/(2*l-1)/(l+1-m)/(l+1+m))*pilm{lp1-1,mp1};
        taulm{lp1+1,mp1}=(l+1)*ct.*pilm{lp1+1,mp1}-(l+1+m)*sqrt((2*(l+1)+1)*(l+1-m)/(2*(l+1)-1)/(l+1+m)).*pilm{lp1,mp1};
    end
end    
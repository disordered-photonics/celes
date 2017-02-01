%======================================================================
%> @brief Map the SVWF coefficients and the sphere number to the
%> multi-index jmult
%>
%> @param       jS (int): Particle number
%> @param       tau (int): SVWF polarization
%> @param       l (int): SVWF degree (polar quantum number)
%> @param       m (int): SVWF order (azimuthal quantum number)
%> @param       lmax (int): SVWF degree cutoff
%> @retval      jmult(int): multi-index
%======================================================================
function jmult = multi2single_index(jS,tau,l,m,lmax)
% jmult = multi2single_index(jS,tau,l,m,lmax)
%
% Return one index subsuming over the indices:
% - jS=1,2,...  (sphere number)
% - tau=1,2     (polarization: 1=TE commonly denoted M, 2=TM commonly denoted N)
% - l=1,..,lmax (polar eigenvalue)
% - m=-l,..,l   (azimuthal eigenvalue)
%
% The multiindex is counted according to the following scheme:
% Count first m, then l, then tau, then jS.
%
%     jmult   |   jS,tau,l,m
%    -----------------------
%         1   |   1,1,1,-1   
%         2   |   1,1,1,0    
%         3   |   1,1,1,1    
%         4   |   1,1,2,-2  
%         5   |   1,1,2,-1   
%         .   |    ....      
%         .   |   1,1,lmax,lmax
%         .   |   1,2,1,-1
%         .   |   1,2,1,0
%         .   |    ....
%         .   |   1,2,lmax,lmax
%         .   |   2,1,1,-1
%         .   |    ....
%   jmult_max |   1,2,1,0

jmult = (jS-1)*2*lmax*(lmax+2)+(tau-1)*lmax*(lmax+2)+(l-1)*(l+1)+m+l+1;
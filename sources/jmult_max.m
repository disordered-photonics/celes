%======================================================================
%> @brief Return the total number of indices
%>
%> @param       nSph (int): Number of particles
%> @param       lmax (int): SVWF degree cutoff
%> @retval      jmax(int): number of indices
%======================================================================

function jmax = jmult_max(nSph,lmax)
% jmax = jmult_max(nSph,lmax)
% What is the dimension of the total index vector?
jmax = 2*nSph*lmax*(lmax+2);
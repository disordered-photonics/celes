%======================================================================
%> @brief Returns the associated Legendre function constructed from a
%> coefficient table. This function is only used for testing.
%>
%> @param ct (float): cos theta
%> @param st (float): sin theta
%> @param l (int): degree of Legendre function
%> @param m (int): order of Legendre function
%> @param coeff_tab (3-dim float array): coefficient table as genereated by
%> Plm_coefficients().
%>
%> @retval Plm (float): Legendre function value
%======================================================================
function Plm = check_legendre(ct,st,l,m,coeff_tab)

Plm = 0;
jj=0;
for lbd=(l-m):-2:0
    jj=jj+1;
    Plm = Plm + st^m * ct^lbd * coeff_tab(l+1,m+1,jj);
end




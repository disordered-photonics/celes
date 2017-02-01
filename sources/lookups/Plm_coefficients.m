%======================================================================
%> @brief Coefficients for the calculation of the associated Legendre
%> functions (requires Matlab's symbolic toolbox)
%>
%> @param lmax (int)
%>
%> @retval Plm_coeff_table 
%> (single prec. float array of size lmax+1 x lmax+1 x ceil(lmax/2))
%> Plm_coeff_table can be used to reconstruct the actual Legendre functions
%> for example by using the function check_legendre()
%======================================================================
function Plm_coeff_table = Plm_coefficients(lmax)

syms ct
syms st
plm = legendre_normalized_trigon(ct,st,lmax);

Plm_coeff_table = zeros(lmax+1,lmax+1,ceil(lmax/2));

for l=0:lmax
    for m=0:l
        [cf,T]=coeffs(plm{l+1,m+1});
        Plm_coeff_table(l+1,m+1,1:length(cf))=single(cf);
    end
end
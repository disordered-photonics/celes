%======================================================================
%> @brief Generate tabulated coefficients for the SVWF translation operator
%>
%> @param lmax (int): truncation degree of the expansion in SVWF
%>
%> @retval translation (struct): structure containing the a5 and b5 
%> coefficients for the SVWF addition theorem, see \ref theory.
%======================================================================
function translation = translation_table_ab(lmax)

jmax = jmult_max(1,lmax);

translation.ab5 = zeros(jmax,jmax,2*lmax+1,'single');

for tau1=1:2
    for l1=1:lmax
        for m1=-l1:l1
            j1=multi2single_index(1,tau1,l1,m1,lmax);
            for tau2=1:2
                for l2=1:lmax
                    for m2=-l2:l2
                        j2=multi2single_index(1,tau2,l2,m2,lmax);
                        for p=0:2*lmax
                            if tau1==tau2
                                translation.ab5(j1,j2,p+1) = (1i)^(abs(m1-m2)-abs(m1)-abs(m2)+l2-l1+p) * (-1)^(m1-m2) ...
                                    * sqrt((2*l1+1)*(2*l2+1)/(2*l1*(l1+1)*l2*(l2+1))) * (l1*(l1+1)+l2*(l2+1)-p*(p+1)) * sqrt(2*p+1) ...
                                    * Wigner3j([l1,l2,p],[m1,-m2,-m1+m2]) * Wigner3j([l1,l2,p],[0,0,0]);
                            elseif p>0
                                translation.ab5(j1,j2,p+1) = (1i)^(abs(m1-m2)-abs(m1)-abs(m2)+l2-l1+p) * (-1)^(m1-m2) ...
                                    * sqrt((2*l1+1)*(2*l2+1)/(2*l1*(l1+1)*l2*(l2+1))) * sqrt((l1+l2+1+p)*(l1+l2+1-p)*(p+l1-l2)*(p-l1+l2)*(2*p+1)) ...
                                    * Wigner3j([l1,l2,p],[m1,-m2,-m1+m2]) * Wigner3j([l1,l2,p-1],[0,0,0]);
                            end
                        end
                    end
                end
            end
        end
    end
end
end
function W = rotation_translation_rotation(k, receiver_position, emitter_position, lmax)

dr = receiver_position - emitter_position;  % in contrast do [Doicu], according to which it would need to be emitter_position - receiver_position
r = norm(dr);
theta = acos(dr(3) / r);
phi = atan2(dr(2), dr(1));

[A, B] = axial_translation_coefficients(lmax, k*r);

D1_lookup = {};
D2_lookup = {};
for m1 = -lmax:lmax
    for m2 = -lmax:lmax
        D1_lookup{m1+lmax+1, m2+lmax+1} = wigner_capital_D(lmax, m1, m2, phi, theta, 0);
        D2_lookup{m1+lmax+1, m2+lmax+1} = wigner_capital_D(lmax, m1, m2, 0, -theta, -phi);
    end
end

W = zeros(jmult_max(1, lmax));

for tau1 = 1:2
    for l1 = 1:lmax
        for m1 = -l1:l1
            j1 = multi2single_index(1, tau1, l1, m1, lmax);
            for tau2 = 1:2
                for l2=1:lmax
                    for m2=-l2:l2
                        j2 = multi2single_index(1, tau2, l2, m2, lmax);
                        for m=-l1:l1
                            D1 = D1_lookup{m1+lmax+1, m+lmax+1}(l1);
                            D2 = D2_lookup{m+lmax+1, m2+lmax+1}(l2);
                            if tau1==tau2
                                AB = A(abs(m)+1, l1, l2);
                            else
                                if m<0
                                    AB = -B(-m+1, l1, l2);
                                else
                                    AB = B(m+1, l1, l2);
                                end
                            end
                            if abs(m)<=min(l1, l2)
                                W(j2, j1) = W(j2, j1) + D1 * AB * D2;  % remember that W is the transpose of A or B
%                                 W(j2, j1) = D2;  % remember that W is the transpose of A or B
                            end
                        end
                    end
                end
            end
        end
    end
end
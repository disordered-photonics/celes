function [A, B, C] = axial_translation_coefficients(lmax, kz)
% returns matrices A and B as defined in [Doicu B.72f (page 285)] 
% For example, A(m+1, l, lprime) contains the  coefficients A^3_{m l, m l'}


lprime_max = 2*lmax+1;

C = NaN([lmax+1, lmax+1, lprime_max+1]);  % indices: m+1, l+1, lprime+1

for m = 0:lmax
    for lprime = 0:lprime_max
        for l = 0:lmax
            if m>l || m>lprime
                C(m+1, l+1, lprime+1) = 0;
            end
        end
    end
end

% [Doicu Appendix B (page 280)]
for lprime = 0:lprime_max
    C(1, 1, lprime+1) = (-1)^lprime * sqrt(2*lprime+1) * sph_bessel(3,lprime,kz);
end

for m = 1:lmax
    for lprime = m:(lprime_max-1)
        C(m+1, m+1, lprime+1) = sqrt((2*m+1)/(2*m)) * ( ...
            sqrt( (lprime + m - 1) * (lprime + m) / (2*lprime-1) / (2*lprime+1)) * C(m, m, lprime) ...
            + sqrt( (lprime-m+1)*(lprime-m+2) / (2*lprime+1) / (2*lprime+3) ) * C(m, m, lprime+2) );
    end
end

% [Doicu B.65 (page 279)]
for m = 0:lmax
    for l = m:(lmax-1)  % we write the coefficient C_{m l+1, m lprime}
        for lprime = m:(lprime_max-(l+1))
            
            if l == 0
                term1 = 0;
            else
                term1 = sqrt( (l-m)*(l+m) / (2*l-1) / (2*l+1) ) * C(m+1, l, lprime+1);
            end
            
            if lprime == 0
                term2 = 0;
            else
                term2 = sqrt( (lprime-m)*(lprime+m) / (2*lprime-1) / (2*lprime+1)) * C(m+1, l+1, lprime);
            end
            
            term3 = - sqrt((lprime-m+1)*(lprime+m+1)/(2*lprime+1)/(2*lprime+3)) * C(m+1, l+1, lprime+2);
            
            C(m+1, l+2, lprime+1) = (term1 + term2 + term3) / sqrt((l-m+1)*(l+m+1)/(2*l+1)/(2*l+3));
            
        end
    end
end

% C = C(:, :, 1:(lmax+2));

A = zeros([lmax+1, lmax, lmax]);  % indices m+1, l, lprime
B = zeros([lmax+1, lmax, lmax]);

% [Doicu Appendix B (page 285)]
for m = 0:lmax
    for lprime = max([1, m]):lmax
        for l = max([1, m]):lmax
            A(m+1, l, lprime) = sqrt((lprime*(lprime+1)) / (l*(l+1))) * ...
                (C(m+1, l+1, lprime+1) ...
                + kz/(lprime+1) * sqrt((lprime-m+1)*(lprime+m+1) / ((2*lprime+1)*(2*lprime+3))) * C(m+1, l+1, lprime+2) ...
                + kz/lprime * sqrt((lprime-m)*(lprime+m) / ((2*lprime+1)*(2*lprime-1))) * C(m+1, l+1, lprime));
            
            B(m+1, l, lprime) = 1i * kz * m/sqrt(l*lprime*(l+1)*(lprime+1)) * C(m+1, l+1, lprime+1);
        end
    end
end

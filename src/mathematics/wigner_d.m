function d = wigner_d(lmax,m,m_prime,beta)
% Returns an array of wigner d symbols d^l_{m,m_prime}(beta) for polar
% angle beta and azimuthal quantum numbers m, m_prime.
% The l-th entry of the result vector corresponds to the polar quantum
% number l.

if beta < 0  % [Doicu, B.41 (page 272)]
    tmp = m;
    m = m_prime;
    m_prime = tmp;
end

smin = max(abs(m),abs(m_prime));  % [Mishchenko, B.22-24 (page 364)]
if m_prime >= m  % [Mishchenko, B.16 (page 364)]
    zeta = 1;
else
    zeta = (-1)^(m-m_prime);
end

x = cos(beta);  % [Mishchenko (page 363)]

if m==0 && m_prime==0
    d = legendre_polynomial(lmax, x);  % [Mishchenko, B.27 (page 365)]
    d = d(2:end);
else
    
    d = zeros(lmax+2,1);  % n=-1 to n=n_max
    
    d(smin+2) = ...    % [Mishchenko, B.24 (page 364)]
        zeta*2^(-smin) ...
        * sqrt(factorial(2*smin) / (factorial(abs(m-m_prime))*factorial(abs(m+m_prime)))) ...
        * (1-x)^(abs(m-m_prime)/2) * (1+x)^(abs(m+m_prime)/2);
    
    for s = smin:(lmax-1)  % [Mishchenko, B.22 (page 365)]
        d(s+3) = ((2*s+1) * (s*(s+1)*x-m*m_prime) * d(s+2) ...
            -(s+1) * sqrt(s^2-m^2) * sqrt(s^2-m_prime^2) * d(s+1)) ...
            / (s * sqrt((s+1)^2-m^2) * sqrt((s+1)^2-m_prime^2));
    end
    
    d = d(3:end); % cut off d=-1,0
end
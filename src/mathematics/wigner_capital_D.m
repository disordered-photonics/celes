function D = wigner_capital_D(lmax, m, m_prime, alpha, beta, gamma)
% Returns an array of wigner D symbols D^l_{m,m_prime}(apha, beta, gamma) 
% for euler angles alpha, beta, gamma and azimuthal quantum numbers m, m_prime.
% The l-th entry of the result vector corresponds to the polar quantum
% number l.

% [Doicu, B.36 (page 271)]
if m >= 0 && m_prime >= 0
    delta = 1;
elseif m >= 0 && m_prime < 0
    delta = (-1)^m_prime;
elseif m < 0 && m_prime >= 0
    delta = (-1)^m;
elseif m < 0 && m_prime < 0
    delta = (-1)^(m+m_prime);
end

d = wigner_d(lmax, m, m_prime, beta);

% [Doicu, B.34 (page 271)]
D = (-1)^(m+m_prime) * delta * exp(1i* (m*alpha + m_prime*gamma)) * d;
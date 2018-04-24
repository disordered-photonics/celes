%  Copyright (c) 2017, Amos Egel (KIT), Lorenzo Pattelli (LENS)
%                      Giacomo Mazzamuto (LENS)
%  All rights reserved.
%
%  Redistribution and use in source and binary forms, with or without
%  modification, are permitted provided that the following conditions are met:
%
%  * Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
%
%  * Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
%
%  * Neither the name of the copyright holder nor the names of its
%    contributors may be used to endorse or promote products derived from
%    this software without specific prior written permission.
%
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
%  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%  POSSIBILITY OF SUCH DAMAGE.

%===============================================================================
%> @brief Evaluate the scattered near field
%>
%> @param   simulation (celes_simulation object): simulation for
%>          which the scattered field is to be computed
%>
%> @retval  E (Nx3 float array): electric field in the format [Ex;Ey;Ez]
%>          Each column correspond to one field point specified in
%>          simulation.output.fieldPoints
%===============================================================================
function E = compute_scattered_field(simulation)

msg='';

% spherical bessel lookup
Rlast = simulation.output.fieldPoints-simulation.input.particles.positionArray(simulation.input.particles.number,:);
rlast = sqrt(sum(Rlast.^2,2));
rmax = max(rlast); % estimate for maximal distance
resol = simulation.numerics.particleDistanceResolution; % resolution of lookup
ri = 0:resol:(rmax+resol);
hi = cell(simulation.numerics.lmax,1);
dhi = cell(simulation.numerics.lmax,1);
for l = 1:simulation.numerics.lmax
    hi{l} = simulation.numerics.deviceArray(single(sph_bessel(3,l,simulation.input.k_medium * ri)));
    dhi{l} = simulation.numerics.deviceArray(single(dx_xz(3,l,simulation.input.k_medium * ri)));
end
ri = simulation.numerics.deviceArray(single(ri));

% initialize field
E = simulation.numerics.deviceArray(zeros(size(simulation.output.fieldPoints),'single'));

for jS = 1:simulation.input.particles.number
    fprintf(1,repmat('\b',1,length(msg)));
    msg = sprintf(' sphere %i of %i',jS,simulation.input.particles.number);
    fprintf(1,msg);

    kM = simulation.input.k_medium;

    % relative positions
    R = simulation.numerics.deviceArray(single(simulation.output.fieldPoints-simulation.input.particles.positionArray(jS,:)));
    r = sqrt(sum(R.^2,2));
    e_r = R./r;
    ct = e_r(:,3);
    ct(ct < -1) = -1; % check for rounding errors
    ct(ct > 1) = 1;
    st = sqrt(1-ct.^2);
    phi = atan2(R(:,2),R(:,1));

    e_theta = [ct.*cos(phi),ct.*sin(phi),-st];
    e_phi = [-sin(phi),cos(phi),zeros(size(phi),'like',phi)];

    kr = kM*r;

    % check if bessel lookup is sufficient
    if max(r) > ri(end)
        rnew = (ri(end)+resol):resol:((max(r)+resol)*1.4);
        for l = 1:simulation.numerics.lmax
            hi{l} = [hi{l},simulation.numerics.deviceArray(single(sph_bessel(3,l,simulation.input.k_medium* gather(rnew))))];
            dhi{l} = [dhi{l},simulation.numerics.deviceArray(single(dx_xz(3,l,simulation.input.k_medium* gather(rnew))))];
        end
        ri = simulation.numerics.deviceArray(single([ri,rnew]));
    end
    
    % spherical functions
    [p_all] = legendre_normalized_trigon(ct,st,simulation.numerics.lmax);
    [pi_all,tau_all] = spherical_functions_trigon(ct,st,simulation.numerics.lmax);

    for l = 1:simulation.numerics.lmax
        z = interp1(ri,hi{l},r,'linear');
        dxxz = interp1(ri,dhi{l},r,'linear');

        for m = -l:l
            P_lm = p_all{l+1,abs(m)+1};
            pi_lm = pi_all{l+1,abs(m)+1};
            tau_lm = tau_all{l+1,abs(m)+1};
            eimphi = exp(1i*m*phi);
            for tau = 1:2
                n = multi2single_index(1,tau,l,m,simulation.numerics.lmax);
                % SVWFs
                if tau == 1  %select M
                    fac1 = 1/sqrt(2*l*(l+1)) * z .* eimphi;
                    fac2 = (1i*m*pi_lm).*e_theta - tau_lm.*e_phi;
                    N = fac1.*fac2;
                else %select N
                    fac1 = 1/sqrt(2*l*(l+1)) * eimphi;
                    term1 = (l*(l+1)*z./kr.*P_lm).*e_r;
                    term22 = tau_lm.*e_theta;
                    term23 = (1i*m*pi_lm).*e_phi;
                    term2 = (dxxz./kr).*(term22+term23);
                    N = fac1.*(term1+term2);
                end
                E = E + simulation.tables.scatteredFieldCoefficients(jS,n) * N;
            end
        end
    end
end
end

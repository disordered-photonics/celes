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
%> @brief Evaluate the internal near field, i.e., the field inside the
%>        scatterers.
%>
%> @param   simulation (celes_simulation object): simulation for
%>          which the scattered field is to be computed
%>
%> @retval  E, H (Nx3 float arrays): electric and magnetic fields in the
%>          format [Ex;Ey;Ez]. Each column correspond to one field point
%>          specified in simulation.output.fieldPoints
%===============================================================================
function [E,H,internal_indices] = compute_internal_field(simulation)

msg='';

% initialize field
E = zeros(size(simulation.output.fieldPoints),'single');
H = zeros(size(simulation.output.fieldPoints),'single');
internal_indices = [];

kM = simulation.input.k_medium;
omega = simulation.input.omega;

for jS=1:simulation.input.particles.number
    fprintf(1,repmat('\b',1,length(msg)));
    msg = sprintf(' sphere %i of %i', jS, simulation.input.particles.number);
    fprintf(1,msg);
    kS = simulation.input.k_particle(jS);
    nu = 1;

    % relative positions
    R = simulation.output.fieldPoints - simulation.input.particles.positionArray(jS,:);
    r2 = sum(R.^2,2);
    intidx = find(r2 < simulation.input.particles.radiusArray(jS)^2);
    Rint = R(intidx,:);

    for l = 1:simulation.numerics.lmax
        for m = -l:l
            [M, N] = SVWF(kS,Rint,nu,l,m); % not using interpolation ...

            n1 = multi2single_index(1,1,l,m,simulation.numerics.lmax);
            n2 = multi2single_index(1,2,l,m,simulation.numerics.lmax);

            b_to_c = T_entry(1,l,kM,kS,simulation.input.particles.radiusArray(jS),'internal')/...
                     T_entry(1,l,kM,kS,simulation.input.particles.radiusArray(jS),'scattered');
            a_to_d = T_entry(2,l,kM,kS,simulation.input.particles.radiusArray(jS),'internal')/...
                     T_entry(2,l,kM,kS,simulation.input.particles.radiusArray(jS),'scattered');

            E(intidx,:) = E(intidx,:) + ...
                simulation.tables.scatteredFieldCoefficients(jS,n1) * b_to_c * M + ...
                simulation.tables.scatteredFieldCoefficients(jS,n2) * a_to_d * N;
            H(intidx,:) = H(intidx,:) -1i*kS/omega *( ...
                simulation.tables.scatteredFieldCoefficients(jS,n1) * b_to_c * N + ...
                simulation.tables.scatteredFieldCoefficients(jS,n2) * a_to_d * M);
        end
    end
    internal_indices = [internal_indices; intidx];
end
end

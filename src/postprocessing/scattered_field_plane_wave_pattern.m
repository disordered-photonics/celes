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
%> @brief Compute the expansion of the scattered field in plane vector wave
%>        functions
%>
%> @param   simulation (celes_simulation object): simulation for
%>          which the scattered field pwp shall be computed
%>
%> @retval  pwp (1x2 cell array): pwp{1} and pwp{2} each contain a
%>          celes_planeWavePattern for the scattered field in TE and
%>          in TM polarization, respectively
%===============================================================================
function pwp = scattered_field_plane_wave_pattern(simulation)

msg='';

lmax = simulation.numerics.lmax;
nmax = simulation.numerics.nmax;
betaArray = simulation.numerics.polarAnglesArray;
cb = cos(betaArray);
sb = sin(betaArray);
alphaArray = simulation.numerics.azimuthalAnglesArray;
alphaDim = length(alphaArray);
betaDim = length(betaArray);

for pol = 2:-1:1
    pwp{pol} = celes_planeWavePattern('polarAngles', betaArray, ...
                                      'azimuthalAngles', alphaArray, ...
                                      'k', simulation.input.k_medium, ...
                                      'expansionCoefficients', simulation.numerics.deviceArray(zeros(alphaDim,betaDim,'single')) ... % dima x dimb
                                      );
end

kxGrid = pwp{1}.kxGrid;
kyGrid = pwp{1}.kyGrid;
kzGrid = pwp{1}.kzGrid;

[pilm,taulm] = spherical_functions_trigon(cb,sb,lmax); % 1 x dimb

for pol = 2:-1:1
    B{pol} = simulation.numerics.deviceArray(zeros(nmax,betaDim,'single')); % dimn x dimb
end

for tau = 1:2
    for l = 1:lmax
        for m = -l:l
            n = multi2single_index(1,tau,l,m,lmax);
            for pol = 1:2
                B{pol}(n,:) = transformation_coefficients(pilm,taulm,tau,l,m,pol);
            end
        end
    end
end

eima = simulation.numerics.deviceArray(zeros(alphaDim,nmax,'single')); % dima x dimn
for tau = 1:2
    for l = 1:lmax
        for m = -l:l
            n = multi2single_index(1,tau,l,m,lmax);
            eima(:,n) = exp(1i*m*alphaArray);
        end
    end
end

for jS = 1:simulation.input.particles.number
    fprintf(1,repmat('\b',1,length(msg)));
    msg = sprintf(' sphere %i of %i',jS,simulation.input.particles.number);
    fprintf(1,msg);

    emnirk = exp(-1i*(simulation.input.particles.positionArray(jS,1)*kxGrid+ ...
                      simulation.input.particles.positionArray(jS,2)*kyGrid+ ...
                      simulation.input.particles.positionArray(jS,3)*kzGrid)); % dima x dimb
    beima = eima.*simulation.tables.scatteredFieldCoefficients(jS,:); % dima x dimn
    for pol = 1:2
        pwp{pol}.expansionCoefficients = pwp{pol}.expansionCoefficients+ ...
                                        (beima * B{pol}).* emnirk / (2*pi);
    end
end
end

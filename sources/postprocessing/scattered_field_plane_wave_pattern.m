%======================================================================
%> @brief Compute the expansion of the scattered field in plane vector wave
%> functions
%>
%> @param       simulation (celes_simulation object): simulation for
%>              which the scattered field pwp shall be computed
%>
%> @retval      pwp (1x2 cell array): pwp{1} and pwp{2} each contain a
%>              celes_planeWavePattern for the scattered field in TE and 
%>              in TM polarization, respectively
%======================================================================
function pwp = scattered_field_plane_wave_pattern(simulation)

msg='';

lmax = simulation.numerics.lmax;
nmax = simulation.numerics.nmax;
betaArray = simulation.numerics.polarAnglesArray;
alphaArray = simulation.numerics.azimuthalAnglesArray;
alphaDim = length(alphaArray);
betaDim = length(betaArray);

for pol=1:2
    pwp{pol} = celes_planeWavePattern(betaArray, alphaArray, simulation.input.k_medium);
end

kxGrid = pwp{1}.kxGrid;
kyGrid = pwp{1}.kyGrid;
kzGrid = pwp{1}.kzGrid;

[pilm,taulm] = spherical_functions_angular(betaArray,lmax);  % 1 x dimb

for pol=1:2
    B{pol} = simulation.numerics.deviceArray(zeros(nmax,betaDim,'single'));  % dimn x dimb
end

for tau=1:2
    for l=1:lmax
        for m=-l:l
            n=multi2single_index(1,tau,l,m,lmax);
            for pol=1:2
                B{pol}(n,:) = transformation_coefficients(pilm,taulm,tau,l,m,pol); 
            end
        end
    end
end

eima = simulation.numerics.deviceArray(zeros(alphaDim,nmax,'single')); % dima x dimn
for tau=1:2
    for l=1:lmax
        for m=-l:l
            n=multi2single_index(1,tau,l,m,lmax);
            eima(:,n) = exp(1i*m*alphaArray);
        end
    end
end

for pol=1:2
    pwp{pol}.expansionCoefficients = simulation.numerics.deviceArray(zeros(alphaDim,betaDim,'single'));  % dima x dimb
end

for jS=1:simulation.input.particles.number
    fprintf(1,repmat('\b',1,length(msg)));
    msg = sprintf(' sphere %i of %i',jS,simulation.input.particles.number);
    fprintf(1,msg);

    emnirk = exp(-1i*(simulation.input.particles.positionArray(jS,1)*kxGrid + simulation.input.particles.positionArray(jS,2)*kyGrid + simulation.input.particles.positionArray(jS,3)*kzGrid));  % dima x dimb
    beima = bsxfun(@times,eima,simulation.tables.scatteredFieldCoefficients(jS,:));  % dima x dimn
    for pol=1:2
        pwp{pol}.expansionCoefficients = pwp{pol}.expansionCoefficients + (beima * B{pol}).* emnirk / (2*pi);
    end
end

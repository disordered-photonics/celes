load('my_simulation_500.mat')

edgeSize = 800;

[partitioning,centerPositions,offsetIdcs] = make_particle_partion(simulation.input.particles.positionArray,[edgeSize,edgeSize,edgeSize]);

lmaxParticle = 3;
lmaxFmm = 10;
ab5_table = translation_table_ab_parallel(lmaxFmm);

fmm = celes_FMM1_structure;

fmm.edgeSize = edgeSize;
fmm.lmax = lmaxFmm;
fmm.centerPositions = centerPositions;
fmm.particleIdcs = partitioning;
fmm.boxOffsetIdcs = offsetIdcs;
fmm.numberOfBoxes = length(partitioning);
fmm = fmm.prepare(simulation.input.particles.positionArray,lmaxParticle,ab5_table,simulation.input.k_medium);

x = simulation.tables.initialFieldCoefficients;

[y3,fmm]=fmm.couplingMatrixMultiply(x);
y4=coupling_matrix_multiply(simulation,x);


%% test aggregation

% actual particle field from box 1:
jbox=1;
E = [0,0,0];
R= [0,0,0];
for jjp = 1:length(fmm.particleIdcs{jbox})
    jp = fmm.particleIdcs{jbox}(jjp);
    for tau=1:2
        for l=1:lmaxParticle
            for m=-l:l
                n=multi2single_index(1,tau,l,m,lmaxParticle);
                b=x(jp,n);
                E = E + b*SVWF(simulation.input.k_medium,R-simulation.input.particles.positionArray(jp,:),3,tau,l,m);
            end
        end
    end
end

%% field from aggregation
E2 = [0,0,0];
for tau=1:2
    for l=1:fmm.lmax
        for m=-l:l
            n=multi2single_index(1,tau,l,m,fmm.lmax);
            b=fmm.boxScatteredFieldCoefficients{jbox}(n);
            E2 = E2 + b*SVWF(simulation.input.k_medium,R-fmm.centerPositions{jbox},3,tau,l,m);
        end
    end
end

%% field from disaggregation:
E4 = [0,0,0];
a = fmm.particleDisaggregationMatrices{jbox} * fmm.boxIncomingFieldCoefficients{jbox};
a = reshape(a,length(fmm.particleIdcs{jbox}),[]);
jjp = 1;%:length(fmm.particleIdcs{jbox})
jp = fmm.particleIdcs{jbox}(jjp);
R = simulation.input.particles.positionArray(jp,:)+[1,1,1]*1e-2;
for tau=1:2
    for l=1:lmaxParticle
        for m=-l:l
            n=multi2single_index(1,tau,l,m,lmaxParticle);
            E4 = E4 + a(jjp,n)*SVWF(simulation.input.k_medium,R-simulation.input.particles.positionArray(jp,:),1,tau,l,m);
        end
    end
end

%% test disaggregation
% actual field from box 1:

E3 = [0,0,0];
for tau=1:2
    for l=1:fmm.lmax
        for m=-l:l
            n=multi2single_index(1,tau,l,m,fmm.lmax);
            a=fmm.boxIncomingFieldCoefficients{jbox}(n);
            E3 = E3 + a*SVWF(simulation.input.k_medium,R-fmm.centerPositions{jbox},1,tau,l,m);
        end
    end
end

%% field via neighborhood coupling matrix

E6 = [0,0,0];
jjnb = 2;
jnb = fmm.neighborBoxIdcs{jbox}(jjnb);

xneigh = x(fmm.particleIdcs{jnb},:);
a = fmm.neighborhoodCouplingMatrices{jbox}{jnb} * xneigh(:);
a = reshape(a,length(fmm.particleIdcs{jbox}),[]);
jjp = 1;%:length(fmm.particleIdcs{jbox})
jp = fmm.particleIdcs{jbox}(jjp);
R = simulation.input.particles.positionArray(jp,:)+[1,1,1]*1e-2;

for tau=1:2
    for l=1:lmaxParticle
        for m=-l:l
            n=multi2single_index(1,tau,l,m,lmaxParticle);
            E6 = E6 + a(jjp,n) * SVWF(simulation.input.k_medium,R-simulation.input.particles.positionArray(jp,:),1,tau,l,m);
        end
    end
end


%% test neighborhood coupling
% actual fieldfrom neighborhood 1


E5 = [0,0,0];

for jjp = 1:length(fmm.particleIdcs{jnb})
    jp = fmm.particleIdcs{jnb}(jjp);
    for tau=1:2
        for l=1:lmaxParticle
            for m=-l:l
                n=multi2single_index(1,tau,l,m,lmaxParticle);
                b=x(jp,n);
                E5 = E5 + b*SVWF(simulation.input.k_medium,R-simulation.input.particles.positionArray(jp,:),3,tau,l,m);
            end
        end
    end
end

%% test box coupling

E7 = [0,0,0];
jbox1 = 1;
jbox2 = 50;
R = fmm.centerPositions{jbox1} - [1,1,1]*1e-2;

for tau=1:2
    for l=1:fmm.lmax
        for m=-l:l
            n=multi2single_index(1,tau,l,m,fmm.lmax);
            b=fmm.boxScatteredFieldCoefficients{jbox2}(n);
            E7 = E7 + b * SVWF(simulation.input.k_medium,R-fmm.centerPositions{jbox2},3,tau,l,m);
        end
    end
end

%%
E8 = [0,0,0];

relativeOffsetIdcs = fmm.boxOffsetIdcs{jbox1} - fmm.boxOffsetIdcs{jbox2};
[~,idx] = ismember(relativeOffsetIdcs,fmm.boxRelativeOffsetIdcs,'rows');
a = fmm.boxCouplingMatrices(:,:,idx) * fmm.boxScatteredFieldCoefficients{jbox2};

for tau=1:2
    for l=1:fmm.lmax
        for m=-l:l
            n=multi2single_index(1,tau,l,m,fmm.lmax);
            E8 = E8 + a(n) * SVWF(simulation.input.k_medium,R-fmm.centerPositions{jbox1},1,tau,l,m);
        end
    end
end









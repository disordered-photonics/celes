load('my_simulation_500.mat')

[partitioning,centerPositions,offsetIdcs] = make_particle_partion(simulation.input.particles.positionArray,[1000,1000,1000]);

lmaxParticle = 3;
lmaxFmm = 5;
ab5_table = translation_table_ab(lmaxFmm);

fmm = celes_FMM1_structure;

fmm.edgeSize = 1000;
fmm.lmax = lmaxFmm;
fmm.centerPositions = centerPositions;
fmm.particleIdcs = partitioning;
fmm.boxOffsetIdcs = offsetIdcs;
fmm.numberOfBoxes = length(partitioning);
fmm = fmm.prepare(simulation.input.particles.positionArray,lmaxParticle,ab5_table.ab5,simulation.input.k_medium)
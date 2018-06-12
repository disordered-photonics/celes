function [Wx] = coupling_matrix_multiply_rtr(simulation,x,varargin)

NS = simulation.input.particles.number;

real_hTab = simulation.h3_table.real_h3;
imag_hTab = simulation.h3_table.imag_h3;

rResol = single(simulation.numerics.particleDistanceResolution) * simulation.input.k_medium;

real_x = gpuArray(single(real(x)));
imag_x = gpuArray(single(imag(x)));

part_pos = gpuArray(single(transpose(simulation.input.particles.positionArray))) * simulation.input.k_medium;

[real_Wx,imag_Wx] = coupling_matrix_multiply_rtr_CUDA(real_x, imag_x, ...
                                                      real_hTab, imag_hTab, ...
                                                      part_pos, int32(NS), rResol);

Wx = real_Wx + 1i*imag_Wx;

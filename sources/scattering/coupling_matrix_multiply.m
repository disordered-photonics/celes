%======================================================================
%> @brief Numerically evaluate the coupling matrix multiply with some
%> coefficient vector. The calculation is done on the GPU
%>
%> This routine is essentially a wrapper for the cuda mex function
%> coupling_matrix_multiply_CUDA
%>
%> @param       simulation (celes_simulation)
%>
%> @param       x (float array): coefficient vector to be multiplied with
%>              the coupling matrix 
%>
%> @retval      Wx (float gpuArray): W*x, where W is the SVWF translation
%>              operator
%======================================================================
function [Wx] = coupling_matrix_multiply(simulation,x)

lmax=simulation.numerics.lmax;
NS = simulation.input.particles.number;

% initialize tables
PlmCoeffTable = gpuArray(single(Plm_coefficients(2*lmax)));

h3_table = spherical_bessel_table(simulation);
real_hTab = gpuArray(single(h3_table.real_h3));
imag_hTab = gpuArray(single(h3_table.imag_h3));
rResol = single(simulation.numerics.particleDistanceResolution);

real_x = gpuArray(single(real(x)));
imag_x = gpuArray(single(imag(x)));

real_ab5Tab=[];
imag_ab5Tab=[];

loopCounter=1;
for tau1=1:2
    for l1=1:lmax
        for m1=-l1:l1
            j1=multi2single_index(1,tau1,l1,m1,lmax);
            for tau2=1:2
                for l2=1:lmax
                    for m2=-l2:l2
                        j2=multi2single_index(1,tau2,l2,m2,lmax);
                        for p=max(abs(m1-m2),abs(l1-l2)+abs(tau1-tau2)):(l1+l2)
                            real_ab5Tab(end+1)=single(real(simulation.tables.translationTable.ab5(j2,j1,p+1)));
                            imag_ab5Tab(end+1)=single(imag(simulation.tables.translationTable.ab5(j2,j1,p+1)));
                            loopCounter=loopCounter+1;
                        end
                    end
                end
            end
        end
    end
end

real_ab5Tab=gpuArray(single(real_ab5Tab));
imag_ab5Tab=gpuArray(single(imag_ab5Tab));
		
part_pos = gpuArray(single(transpose(simulation.input.particles.positionArray)));

% tic
[real_Wx,imag_Wx] = coupling_matrix_multiply_CUDA(real_x,imag_x,real_hTab,imag_hTab,PlmCoeffTable,real_ab5Tab,imag_ab5Tab,part_pos,int32(NS),rResol);
% wait(gpuDevice)
% toc

% tic
% [real_Wx,imag_Wx] = coupling_matrix_multiply_copy_of_CUDA_lmax3(gather(real_x),gather(imag_x),gather(real_hTab),gather(imag_hTab),gather(PlmCoeffTable),gather(real_ab5Tab),gather(imag_ab5Tab),gather(part_pos),gather(int32(NS)),gather(rResol));
% toc
% real_Wx(1)
% imag_Wx(1)

Wx = real_Wx + 1i*imag_Wx;
% norm(Wx)


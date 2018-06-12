function [Wx] = coupling_matrix_multiply_rtr(simulation,x,varargin)

lmax = simulation.numerics.lmax;
NS = simulation.input.particles.number;

Wx = zeros(size(x));

for js1 = 1:NS
    for js2 = 1:NS
        if js1 ~= js2
            W = rotation_translation_rotation(...
                simulation.input.k_medium, ...
                simulation.input.particles.positionArray(js1, :), ...
                simulation.input.particles.positionArray(js2, :), ...
                lmax);
            
            Wx(js1, :) = Wx(js1, :) + transpose( W * transpose(x(js2, :)));
        end
    end
end

%======================================================================
%> @brief Evaluate the internal near field, i.e., the field inside the
%> scatterers.
%>
%> @param       simulation (celes_simulation object): simulation for
%>              which the scattered field is to be computed
%>
%> @retval      E (Nx3 float array): electric field in the format [Ex;Ey;Ez]
%>              Each column correspond to one field point specified in 
%>              simulation.output.fieldPoints
%======================================================================
function [E,internal_indices] = compute_internal_field(simulation)

msg='';

% initialize field
E = zeros(size(simulation.output.fieldPoints),'single');
internal_indices = [];

kM = simulation.input.k_medium;

for jS=1:simulation.input.particles.number
    fprintf(1,repmat('\b',1,length(msg)));
    msg = sprintf(' sphere %i of %i',jS,simulation.input.particles.number);
    fprintf(1,msg);
    
    kS = simulation.input.k_particle;
    nu = 1;
    
    % relative positions
    R = bsxfun(@plus,simulation.output.fieldPoints,-simulation.input.particles.positionArray(jS,:));
    r = sqrt(R(:,1).^2+R(:,2).^2+R(:,3).^2);
    Rint = R(r<simulation.input.particles.radius,:);
    
    for l=1:simulation.numerics.lmax
        for m=-l:l
            for tau=1:2
                n = multi2single_index(1,tau,l,m,simulation.numerics.lmax);
                N = SVWF(kS,Rint,nu,tau,l,m);
                b_to_c = T_entry(tau,l,kM,kS,simulation.input.particles.radius,'internal')/T_entry(tau,l,kM,kS,simulation.input.particles.radius,'scattered');
                E(r<simulation.input.particles.radius,:) = E(r<simulation.input.particles.radius,:) + simulation.tables.scatteredFieldCoefficients(jS,n) * b_to_c * N;
            end
        end
    end
    internal_indices = [internal_indices;find(r<simulation.input.particles.radius)];
end
%======================================================================
%> @brief Rearrange the particle order such each partition is a block in
%> the positionArray
%>
%> @param       simulation (celes_simulation)
%>
%> @param       partitioning (cell array): partitioning according to which
%>              the position array shall be rearranged
%>
%> @retval      simulation (celes_simulation): simulation object with
%>              updated particle order in 
%>              simulation.input.particles.positionArray
%======================================================================
function simulation = sort_particles_by_partition(simulation,partitioning)

partitionArray = cell2mat(partitioning');
simulation.input.particles.positionArray = simulation.input.particles.positionArray(partitionArray,:);

switch simulation.input.particles.disperse
    case 'mono'
    otherwise
        error('only mono disperse supported')
end

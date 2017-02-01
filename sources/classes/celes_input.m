%> @file celes_input.m
% ======================================================================
%> @brief Parameters describing the model to be simulated, i.e., the particles and the excitation
% ======================================================================

classdef celes_input
    
    properties
        %> vacuum wavelength
        wavelength
        
        %> refractive index of the surrounding medium
        mediumRefractiveIndex
        
        %> celes_particles object which contains the parameters that 
        %> specify the particles sizes, positions and refractive indices
        particles
        
        %> celes_initialField object which contains the parameters that 
        %> specify the initial excitation (Gaussian beam)
        initialField
    end
    
    properties (Dependent)
        %> angular frequency (we use c=1)
        omega
        
        %> wavenumber in surrounding medium
        k_medium
        
        %> wavenumber inside particles
        k_particle
    end
    
    methods
        % ======================================================================
        %> @brief Get method for omega
        % ======================================================================
        function value=get.omega(obj)
            value = single(2*pi/obj.wavelength);
        end
        
        % ======================================================================
        %> @brief Get method for k_medium
        % ======================================================================
        function value=get.k_medium(obj)
            value = obj.omega*obj.mediumRefractiveIndex;
        end
        
        % ======================================================================
        %> @brief Get method for k_particle
        % ======================================================================
        function value=get.k_particle(obj)
            switch obj.particles.disperse
                case 'mono'
                    value = obj.omega*obj.particles.refractiveIndex;
                otherwise
                    error( 'action undefined' ) % to be implemented
            end
        end
        
        
    end
end


%> @file celes_numerics.m
% ======================================================================
%> @brief Parameters describing the numerical settings used in the simulation
% ======================================================================

classdef celes_numerics
    properties
        %> maximal polar multipole order
        lmax
        
        %> array of polar angles (between 0 and pi) that specifies the
        %> discretization of angles in the plane wave patterns
        polarAnglesArray
        
        %> array of azimuthal angles (between 0 and 2*pi) that specifies 
        %> the discretization of angles in the plane wave patterns
        azimuthalAnglesArray
        
        %> shall the GPU be used in computations of, e.g., the
        %> preconditioner?
        %> NOTE: In the current version, the GPU is alway used for the
        %> matrix vector products. Future versions might offer the
        %> opportunity to run the whole calculation on the CPU
        gpuFlag = true
        
        %> celes_solver object that contains all the information about
        %> how the linear system is to be solved
        solver = celes_solver
        
        %> resolution of the radial particle distance in the lookup tables
        %> of the translation operator
        particleDistanceResolution = single(10);
    end
    
    properties (Dependent)
        %> number of unknowns per particle
        nmax
    end
    
    methods
        % ======================================================================
        %> @brief Get method for nmax
        % ======================================================================
        function value=get.nmax(obj)
            value=jmult_max(1,obj.lmax);
        end
        
        % ======================================================================
        %> @brief Returns a gpuArray if gpuFlag, and a usual array
        %> otherwise
        %>
        %> @param Array x
        %> @return gpuArray(x) if gpuFlag, otherwise x
        % ======================================================================
        function arr_out = deviceArray(obj,arr_in)
            if obj.gpuFlag
                arr_out=gpuArray(arr_in);
            else
                arr_out=arr_in;
            end
        end
        
    end
end


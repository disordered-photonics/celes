%> @file celes_planeWavePattern.m
% ======================================================================
%> @brief Expansion of a field in plane vector wave functions
% ======================================================================

classdef celes_planeWavePattern

    properties
        %> array of polar angle values (in rad)
        polarAngles        
        
        %> array of azimuthal angle values (in rad)
        azimuthalAngles     
        
        %> wavenumber
        k
        
        %> coefficients of the PVWF expansion
        expansionCoefficients
    end
    
    properties (Dependent)
        %> meshigrid of polar angles
        polarAnglesGrid
        
        %> meshigrid of azimuthal angles
        azimuthalAnglesGrid
        
        %> array of cylindrical radial components of wave vector
        kparArray
        
        %> array of cylindrical z-components of wave vector
        kzArray
        
        %> meshgrid of cylindrical radial components of wave vector
        kparGrid
        
        %> meshgrid of z-components of wave vector
        kzGrid
        
        %> meshgrid of x-components of wave vector
        kxGrid
        
        %> meshgrid of y-components of wave vector
        kyGrid
    end        
    
    methods
        % ======================================================================
        %> @brief Constructor method
        % ======================================================================
        function obj = celes_planeWavePattern(varargin)
            if ~isempty(varargin)
                obj.polarAngles = varargin{1};
                obj.azimuthalAngles = varargin{2};
                obj.k = varargin{3};
            end
        end
        
        % ======================================================================
        %> @brief Get method for polarAnglesGrid
        % ======================================================================
        function value = get.polarAnglesGrid(obj)
            [value,~]=meshgrid(obj.polarAngles,obj.azimuthalAngles);
        end
        
        % ======================================================================
        %> @brief Get method for azimuthalAnglesGrid
        % ======================================================================
        function value = get.azimuthalAnglesGrid(obj)
            [~,value]=meshgrid(obj.polarAngles,obj.azimuthalAngles);
        end

        % ======================================================================
        %> @brief Get method for kparArray
        % ======================================================================
        function value = get.kparArray(obj)
            value = obj.k*sin(obj.polarAngles);
        end
        
        % ======================================================================
        %> @brief Get method for kzArray
        % ======================================================================
        function value = get.kzArray(obj)
            value = obj.k*cos(obj.polarAngles);
        end
        
        % ======================================================================
        %> @brief Get method for kxGrid
        % ======================================================================
        function value = get.kxGrid(obj)
            value = obj.k*sin(obj.polarAnglesGrid).*cos(obj.azimuthalAnglesGrid);
        end
        
        % ======================================================================
        %> @brief Get method for kyGrid
        % ======================================================================
        function value = get.kyGrid(obj)
            value = obj.k*sin(obj.polarAnglesGrid).*sin(obj.azimuthalAnglesGrid);
        end
        
        % ======================================================================
        %> @brief Get method for kparGrid
        % ======================================================================
        function value = get.kparGrid(obj)
            value = obj.k*sin(obj.polarAnglesGrid);
        end
        
        % ======================================================================
        %> @brief Get method for kzGrid
        % ======================================================================
        function value = get.kzGrid(obj)
            value = obj.k*cos(obj.polarAnglesGrid);
        end
        
        % ======================================================================
        %> @brief Add this plane wave pttern to another one
        %>
        %> @param celes_planeWavePattern pwp2
        %> @return celes_planeWavePattern as the sum of this one with pwp2
        % ======================================================================
        function obj = addTo(obj,obj2)
            if any(size(obj.polarAngles)-size(obj2.polarAngles)) || any(size(obj.azimuthalAngles)-size(obj2.azimuthalAngles)) || ~(obj.k==obj2.k)
                error('plane wave patterns are not consistent')
            else
                obj.expansionCoefficients = obj.expansionCoefficients + obj2.expansionCoefficients;
            end
        end
    end
end
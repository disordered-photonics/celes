%  Copyright (c) 2017, Amos Egel (KIT), Lorenzo Pattelli (LENS)
%                      Giacomo Mazzamuto (LENS)
%  All rights reserved.
%
%  Redistribution and use in source and binary forms, with or without
%  modification, are permitted provided that the following conditions are met:
%
%  * Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
%
%  * Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
%
%  * Neither the name of the copyright holder nor the names of its
%    contributors may be used to endorse or promote products derived from
%    this software without specific prior written permission.
%
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
%  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%  POSSIBILITY OF SUCH DAMAGE.

%> @file celes_planeWavePattern.m
% ==============================================================================
%> @brief Expansion of a field in plane vector wave functions
% ==============================================================================

classdef celes_planeWavePattern < matlab.System

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

    properties (SetAccess=private, Hidden)
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
            setProperties(obj,nargin,varargin{:});
            validatePropertiesImpl(obj);
            setupImpl(obj);
        end

        % ======================================================================
        %> @brief Add this plane wave pattern to another one
        %>
        %> @param celes_planeWavePattern pwp2
        %> @return celes_planeWavePattern as the sum of this one with pwp2
        % ======================================================================
        function obj = addTo(obj,obj2)
            if any(size(obj.polarAngles)-size(obj2.polarAngles)) || ...
               any(size(obj.azimuthalAngles)-size(obj2.azimuthalAngles)) || ...
               ~(obj.k==obj2.k)
                error('plane wave patterns are not consistent')
            else
                obj.expansionCoefficients = obj.expansionCoefficients+ ...
                                            obj2.expansionCoefficients;
            end
        end
    end

    methods(Access = protected)
        % ======================================================================
        %> @brief Class implementation
        % ======================================================================
        function setupImpl(obj)
            generateGrids(obj);
        end

        % ======================================================================
        %> @brief Validate class properties
        % ======================================================================
        function validatePropertiesImpl(obj)
            try validateattributes(obj.polarAngles,{'numeric'},{'real','nonnan','finite'})
            catch e, error('invalid polarAngles array: %s', e.message); end
            try validateattributes(obj.azimuthalAngles,{'numeric'},{'real','nonnan','finite'})
            catch e, error('invalid azimuthalAngles array: %s', e.message); end
            try validateattributes(obj.k,{'numeric'},{'nonnan','real','finite','scalar'})
            catch e, error('invalid k value: %s', e.message); end
            try validateattributes(obj.expansionCoefficients,{'gpuArray','numeric'},{'nonnan','finite','2d'})
            catch e, error('invalid expansionCoefficients: %s', e.message); end
        end

        % ======================================================================
        %> @brief Create grids
        % ======================================================================
        function generateGrids(obj)
            obj.kparArray = obj.k*sin(obj.polarAngles);
            obj.kzArray = obj.k*cos(obj.polarAngles);
            [obj.polarAnglesGrid, obj.azimuthalAnglesGrid] = ...
                                meshgrid(obj.polarAngles,obj.azimuthalAngles);

            obj.kparGrid = obj.k*sin(obj.polarAnglesGrid);

            obj.kxGrid = obj.kparGrid.*cos(obj.azimuthalAnglesGrid);
            obj.kyGrid = obj.kparGrid.*sin(obj.azimuthalAnglesGrid);
            obj.kzGrid = obj.k*cos(obj.polarAnglesGrid);
        end

    end
end

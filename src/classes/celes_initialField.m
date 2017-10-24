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

%> @file celes_initialField.m
% ==============================================================================
%> @brief Parameters describing the initial field, e.g. an incoming Gaussian
%>        beam.
% ==============================================================================

classdef celes_initialField < matlab.System

    properties
        %> specify the type of the field:
        %> currently, only 'gaussian beam' is possible
        type  = 'gaussian beam'

        %> amplitude of initial beam
        amplitude = single(1)

        %> incident angle (polar, in radians)
        %> currently, only normal incidence (0) is possible
        polarAngle = single(0)

        %> incident angle (azimuth, in radians)
        azimuthalAngle = single(0)

        %> incident beam polarization ('TE' or 'TM')
        polarization = 'TE'

        %> incident beam width at focal point
        beamWidth

        %> focus of incident beam, [x,y,z]
        focalPoint
    end

    properties (SetAccess=private, Hidden)
        %> incident beam polarization (1 for 'TE', 2 for 'TM')
        pol

        %> is the beam coming at normal incidence,
        %> i.e. is sin(polarAngle)=0?
        normalIncidence
    end

    methods
        % ======================================================================
        %> @brief Class constructor
        % ======================================================================
        function obj = celes_initialField(varargin)
            setProperties(obj,nargin,varargin{:});
            validatePropertiesImpl(obj);
            setupImpl(obj);
        end

        % ======================================================================
        %> @brief Set incident beam polarization (1 for 'TE', 2 for 'TM')
        % ======================================================================
        function setPolIdx(obj)
            switch lower(obj.polarization)
                case 'te'
                    obj.pol = 1;
                case 'tm'
                    obj.pol = 2;
                otherwise
                    error('invalid polarization')
            end
        end

        % ======================================================================
        %> @brief Set normalIncidenceFlag
        % ======================================================================
        function setNormalIncidence(obj)
            obj.normalIncidence = (abs(sin(obj.polarAngle))<1e-5);
        end

    end

    methods(Access = protected)
        % ======================================================================
        %> @brief Class implementation
        % ======================================================================
        function setupImpl(obj)
            setPolIdx(obj)
            setNormalIncidence(obj)
        end

        % ======================================================================
        %> @brief Validate class properties
        % ======================================================================
        function validatePropertiesImpl(obj)
            try validateattributes(obj.type,{'char'},{'nonempty'})
            catch e, error('invalid initialField type: %s', e.message); end
            try validateattributes(obj.amplitude,{'numeric'},{'real','nonnan','finite','scalar'})
            catch e, error('invalid amplitude value: %s', e.message); end
            if strcmpi(obj.type,'gaussian beam')
                try validateattributes(obj.polarAngle,{'numeric'},{'real','nonnan','finite','scalar'})
                catch e, error('invalid polarAngle value: %s', e.message); end
                try validateattributes(obj.azimuthalAngle,{'numeric'},{'real','nonnan','finite','scalar'})
                catch e, error('invalid azimuthalAngle value: %s', e.message); end
                try validateattributes(obj.polarization,{'char'},{'nonempty'})
                catch e, error('invalid polarization type: %s', e.message); end
                try validateattributes(obj.beamWidth,{'numeric'},{'real','nonnan','scalar'})
                catch e, error('invalid beamWidth value: %s', e.message); end
                try validateattributes(obj.focalPoint,{'numeric'},{'real','nonnan','finite','row','numel',3})
                catch e, error('invalid focalPoint array: %s', e.message); end
            else
                error('invalid initialField type')
            end
        end
    end
end

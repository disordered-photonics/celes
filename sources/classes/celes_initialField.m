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
% ======================================================================
%> @brief Parameters describing the initial field, e.g. an incoming Gaussian beam.
% ======================================================================

classdef celes_initialField

    properties
        %> specify the type of the field: 
        %> currently, only 'gaussian beam' is possible
        type  = 'gaussian beam'
        
        %> amplitude of initial beam
        amplitude = single(1)
        
        %> incident angle (polar, in radians)
        %> currently, only normal incidence (0) is possible
        polarAngle  = single(0)
        
        %> incident angle (azimuth, in radians)
        azimuthalAngle = single(0)
        
        %> incident beam polarization ('TE' or 'TM')
        polarization = 'TE'
        
        %> incident beam width at focal point
        beamWidth
        
        %> focus of incident beam, [x,y,z]
        focalPoint
    end
    
    properties(Dependent)
        %> incident beam polarization (1 for 'TE', 2 for 'TM')
        pol
        
        %> is the beam coming at normal incidence,
        %> i.e. is sin(polarAngle)=0?
        normalIncidence
    end
    
    methods
        % ======================================================================
        %> @brief Set method for type
        % ======================================================================
        function inF = set.type(inF,value)
            switch(value)
                case 'gaussian beam'
                    inF.type='gaussian beam';
                otherwise
                    error('this is not a legal initial field type')
            end
        end
        
        % ======================================================================
        %> @brief Set method for polarAngle
        % ======================================================================
        function inF = set.polarAngle(inF,value)
            inF.polarAngle=single(value);
        end
        
        % ======================================================================
        %> @brief Set method for azimthalAngle
        % ======================================================================
        function inF = set.azimuthalAngle(inF,value)
            inF.azimuthalAngle=single(value);
        end
        
        % ======================================================================
        %> @brief Set method for polarization
        % ======================================================================
        function inF = set.polarization(inF,value)
            switch(value)
                case 'TE'
                    inF.polarization='TE';
                case 'TM'
                    inF.polarization='TM';
                otherwise
                    error('this is not a legal initial field polarization')
            end
        end
          
        % ======================================================================
        %> @brief Set method for beamWidth
        % ======================================================================
        function inF = set.beamWidth(inF,value)
            inF.beamWidth=single(value);
        end
        
        % ======================================================================
        %> @brief Set method for focalPoint
        % ======================================================================
        function inF = set.focalPoint(inF,value)
            if length(value)==3
                inF.focalPoint=single(value);
            else
                error('illegal focal point')
            end
        end
                
        % ======================================================================
        %> @brief Get method for pol
        % ======================================================================
        function value = get.pol(inF)
            switch inF.polarization
                case 'TE'
                    value=1;
                case 'TM'
                    value=2;
            end
        end
        
        % ======================================================================
        %> @brief Get method for normalIncidence
        % ======================================================================
        function value = get.normalIncidence(obj)
            value = (abs(sin(obj.polarAngle))<1e-5);
        end
    end
end


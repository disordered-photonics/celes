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


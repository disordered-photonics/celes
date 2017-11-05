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

%> @file celes_output.m
% ==============================================================================
%> @brief Holds the output of a celes_simulation
% ==============================================================================

classdef celes_output < matlab.System

    properties
        %> power flux of the initial field
        initialFieldPower

        %> celes_planeWavePattern object for the scattered field
        scatteredFieldPlaneWavePattern

        %> celes_planeWavePattern object for the total field
        totalFieldPlaneWavePattern

        %> power flux of the total field in +z direction
        totalFieldForwardPower

        %> power flux of the total field in -z direction
        totalFieldBackwardPower

        %> electric (near) field values of the initial field at the points
        %> specified in fieldPoints
        initialField

        %> if true, the paraxial approximation is used for the electric
        %> field evaluation of the initial field, otherwise compute initial
        %> field from angular spectrum representation
        initialFieldUseParaxialApproximation = true

        %> electric (near) field values of the scattered field at the points
        %> specified in fieldPoints
        scatteredField

        %> electric (near) field values of the internal field at the points
        %> specified in fieldPoints
        internalField

        %> indices that correspond to field points that are inside some
        %> sphere
        internalIndices

        %> points at which the electric field is to be evaluated, in the
        %> format [x(:),y(:),z(:)]
        fieldPoints

        %> if the field points lie on a plane (in order to display the
        %> field as an image), fieldPointsArrayDims stores the size of the
        %> image - which can be reconstructed by using
        %> reshape(...,fieldPointsArrayDims)
        fieldPointsArrayDims

        %> array of the relative error of the solution as a function of
        %> iteration step number. for BiCGSTAB algorithm, the half-steps
        %> are counted
        convergenceHistory

        %> time used for the simulation
        totalTime

        %> time used by the solver of the linear system
        solverTime

        %> time used for the preparation of the preconditioner
        preconiditionerPreparationTime

        %> time used for the simulation.run command
        runningTime

        %> time used for the power evaluation
        powerEvaluationTime

        %> time used for the field evaluation
        fieldEvaluationTime
    end

    properties (Dependent)
        %> electric (near) field values of the total field at the points
        %> specified in fieldPoints
        totalField
    end

    methods
        % ======================================================================
        %> @brief Class constructor
        % ======================================================================
        function obj = celes_output(varargin)
            if nargin
                setProperties(obj,nargin,varargin{:});
                validatePropertiesImpl(obj);
            end
        end

        % ======================================================================
        %> @brief Get method for totalField
        % ======================================================================
        function E = get.totalField(obj)
            try
                E = obj.initialField + obj.scatteredField;
                E(obj.internalIndices,:) = ...
                                       obj.internalField(obj.internalIndices,:);
            catch
                E = [];
            end
        end
    end

    methods(Access = protected)
        % ======================================================================
        %> @brief Validate class properties
        % ======================================================================
        function validatePropertiesImpl(obj)
            try validateattributes(obj.fieldPoints,{'numeric'},{'real','nonnan','finite','2d','ncols',3})
            catch e, error('invalid fieldPoints array: %s', e.message); end
            try validateattributes(obj.fieldPointsArrayDims,{'numeric'},{'integer','finite','nonnan','numel',2})
            catch e, error('invalid fieldPointsArrayDims values: %s', e.message); end
        end
    end
end

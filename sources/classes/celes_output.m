%> @file celes_output.m
% ======================================================================
%> @brief Holds the output of a celes_simulation
% ======================================================================

classdef celes_output

    properties
        %> power flux of the initial field
        initialFieldPower
        
        %> celes_planeWavePattern object for the scattered field
        scatteredFieldPlaneWavePattern
        
        %> celes_planeWavePattern object for the total field
        totalFieldPlaneWavePattern
        
        %> power flux of the toal field in +z direction
        totalFieldForwardPower
        
        %> power flux of the toal field in -z direction
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
        %> image - which can be recustructed by using
        %> reshape(...,fieldPointsArrayDims)
        fieldPointsArrayDims
        
        %> array of the relative error of the solution as a function of
        %> iteration step number. for bicgstab algorithm, the half-steps
        %> are counted
        convergenceHistory
        
        %> time used for the simulation
        totalTime
        
        %> time used by the solver of the linear system
        solverTime
        
        %> time used for the preparation of the preconditioner
        preconiditionerPreparationTime
    end
    
    properties (Dependent)
        %> electric (near) field values of the total field at the points 
        %> specified in fieldPoints
        totalField
    end
    
    methods
        % ======================================================================
        %> @brief Get method for totalField
        % ======================================================================
        function E = get.totalField(obj)
            try
                E = obj.initialField + obj.scatteredField;
                E(obj.internalIndices,:) = obj.internalField(obj.internalIndices,:);
            catch
                E = [];
            end
        end
    end
end


%> @file celes_tables.m
% ======================================================================
%> @brief Objects of this class hold large tables of interim results
% ======================================================================

classdef celes_tables

    properties
        %> a table of coefficients needed for the SVWF translation
        translationTable
        
        %> T-matrix of the spheres
        mieCoefficients
        
        %> coefficients of the regular SVWF expansion of the initial
        %> excitation 
        initialFieldCoefficients
        
        %> coefficients of the outgoing SVWF expansion of the scattered
        %> field 
        scatteredFieldCoefficients
    end
    
    properties (Dependent)
        %> right hand side T*aI of linear system M*b=T*aI
        rightHandSide
    end
    
    methods
        % ======================================================================
        %> @brief Get method for rightHandSide
        % ======================================================================
        function TaI = get.rightHandSide(obj)
            TaI = bsxfun(@times,obj.initialFieldCoefficients,obj.mieCoefficients);
        end
    end
end


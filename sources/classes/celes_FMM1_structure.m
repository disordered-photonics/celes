classdef celes_FMM1_structure
    
    properties
        edgeSize
        numberOfBoxes
        centerPositions
        boxOffsetIdcs
        
        particleIdcs
        neighborBoxIdcs = {}
        
        
        lmax
        
        
        particleAggregationMatrices = {}
        particleDisaggregationMatrices = {}
        
        boxRelativeOffsetIdcs = []
        boxCouplingMatrices = {}
        
        neighborhoodCouplingMatrices = {}
        
        boxScatteredFieldCoefficients = {};
        boxIncomingFieldCoefficients = {};
    end
    
    
    methods
        
        function obj = prepare(obj,positionArray,particleLmax,ab5_table,k_medium)
            obj.neighborBoxIdcs = cell(obj.numberOfBoxes,1);
            for jbox1 = 1:obj.numberOfBoxes % receiving box
                for jbox2 = 1:obj.numberOfBoxes % emitting box
                    relativeOffsetIdcs = obj.boxOffsetIdcs{jbox1} - obj.boxOffsetIdcs{jbox2}; % relative offset between receiving and emitting box
                    if max(abs(relativeOffsetIdcs))<=1 &&  sum(abs(relativeOffsetIdcs))<=2 % are boxes neighbors?
                        obj.neighborBoxIdcs{jbox1}(end+1)=jbox2;
                    else
                        obj.boxRelativeOffsetIdcs(end+1,:) = relativeOffsetIdcs;
                    end
                end
            end
            obj.boxRelativeOffsetIdcs = unique(obj.boxRelativeOffsetIdcs,'rows');
            
            obj = assemble_box_coupling_matrices(obj,ab5_table,k_medium);
            obj = assemble_aggregation_matrices(obj,positionArray,particleLmax,ab5_table,k_medium);
            obj = assemble_disaggregation_matrices(obj,positionArray,particleLmax,ab5_table,k_medium);
            obj = assemble_neighborhood_coupling_matrices(obj,ab5_table,positionArray,particleLmax,k_medium);            
        end
        
        function y = couplingMatrixMultiply(obj,x)
            % return y = W*x
            
            y = x-x;
            
            % aggregation
            for jbox = 1:obj.numberOfBoxes
                xBox = x(obj.particleIdcs{jbox},:);
                xBox = xBox(:);
                obj.boxScatteredFieldCoefficients{jbox} = obj.particleAggregationMatrices{jbox} * xBox;
            end
            
            % box coupling
            for jbox1 = 1:obj.numberOfBoxes %receiving box
                obj.boxIncomingFieldCoefficients{jbox1} = zeros(jmult_max(1,obj.lmax),1,'single');
                for jbox2 = 1:obj.numberOfBoxes %emitting box
                    if ~any(jbox2==obj.neighborBoxIdcs{jbox1})
                        relativeOffsetIdcs = obj.boxOffsetIdcs{jbox1} - obj.boxOffsetIdcs{jbox2};
                        [~,idx] = ismember(relativeOffsetIdcs,obj.boxRelativeOffsetIdcs,'rows');
                        obj.boxIncomingFieldCoefficients{jbox1} = obj.boxIncomingFieldCoefficients{jbox1} + obj.boxCouplingMatrices(:,:,idx)*obj.boxScatteredFieldCoefficients{jbox2};
                    end
                end
            end
            
            % disaggregation and neighborhood coupling
            for jbox = 1:obj.numberOfBoxes
                yBox = obj.particleDisaggregationMatrices{jbox} * obj.boxIncomingFieldCoefficients{jbox};
                for jbox2 = obj.neighborBoxIdcs{jbox}
                    xBoxNeighb = x(obj.particleIdcs{jbox2},:);
                    xBoxNeighb = xBoxNeighb(:);
                    yBox = yBox + obj.neighborhoodCouplingMatrices{jbox}{jbox2} * xBoxNeighb;
                end
                y(obj.particleIdcs{jbox},:) = reshape(yBox,length(obj.particleIdcs{jbox}),length(x(1,:)));
            end
        end
    end
end



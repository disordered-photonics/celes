%> @file celes_FMM1_structure.m
% ======================================================================
%> @brief Class which manages the fast multipole algorithm.
% ======================================================================

classdef celes_FMM1_structure
    
    properties
        %> Size of the cubic boxes. The field of particles in one box is
        %> expressed in one common multipole expansion
        edgeSize
        
        %> Total number of boxes
        numberOfBoxes
        
        %> Cell array of positions of the box centers in the format 
        %> {[x1,y1,z1],[x2,y2,z2],...}
        centerPositions
        
        %> Cell array of integer 3-vectors that specify the relative
        %> position of each box with respect to the first box in units of
        %> edgeSize. I.e., 
        %> centerPositions{j}-centerPositions{1}=boxOffsetIdcs{j}*edgeSize
        boxOffsetIdcs
                
        %> Cell array of indices of particles that are located in the
        %> respective box, i.e., if box jb contains particles jp1,jp2,...
        %> then particleIdcs{jb}=[jp1,jp2,...]
        particleIdcs
        
        %> Cell array of indices of boxes that are neighbors of the
        %> respective box.
        neighborBoxIdcs = {}
        
        %> Multipole order truncation for the FMM expansion (should
        %> correspond to the box size)
        lmax
        
        %> Cell array of matrices that map the scattered field coefficients
        %> of the particles inside a box to the scattered field
        %> coefficients of the box
        particleAggregationMatrices = {}
        
        %> Cell array of matrices that map the incoming field coefficients
        %> of the box to the incoming field coefficients of the particles
        particleDisaggregationMatrices = {}
        
        %> Displacement indices [dx1,dy1,dy2;dx2,dy2,dz2;...]
        %> Each row of this array contains integer numbers that correspond
        %> to the difference in box offest indices:
        %> centerPositions{j1}-centerPositions{j2} = boxOffsetIdcs(jj,:)*edgeSize
        %> for some jj.        
        boxRelativeOffsetIdcs = []
        
        %> 3D-Array. boxCouplingMatrices(:,:,idx) corresponds to the matrix
        %> that maps the scattered field coefficients of box j1 to the
        %> incoming field coefficients of box j2 if 
        %> boxOffsetIdcs{j1}-boxOffsetIdcs{j2}=boxRelativeOffsetIdcs(idx,:)
        boxCouplingMatrices = []
        
        %> Cell array of cell arrays of matrices that manage the particle
        %> to particle coupling for neighboring boxes (including the
        %> selected box itself).
        neighborhoodCouplingMatrices = {}
        
        %> Cell array of coefficient vectors
        %> Coefficients of the expansion of the scattered field of the
        %> particles inside the box in terms of multipoles relative to the
        %> box center.
        boxScatteredFieldCoefficients = {};
        
        %> Cell array of coefficient vectors
        %> Coefficients of the expansion of the incoming field of the
        %> particles inside the box in terms of multipoles relative to the
        %> box center. Only the incoming field due to the scattered field
        %> from particles in non-neighboring boxes is regarded.
        boxIncomingFieldCoefficients = {};
    end
    
    
    methods
        % ======================================================================
        %> @brief Prepare 1st level fast multipole algorithm
        %>
        %> @param Optional: b0, initial guess for scattered field
        %> coefficients
        %> @return celes_simulation object with updated initialFieldPower
        % ======================================================================
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
            
            disp('assemble box coupling')
            obj = assemble_box_coupling_matrices(obj,ab5_table,k_medium);
            disp('assemble aggregation')
            obj = assemble_aggregation_matrices(obj,positionArray,particleLmax,ab5_table,k_medium);
            disp('assemble disaggregation')
            obj = assemble_disaggregation_matrices(obj,positionArray,particleLmax,ab5_table,k_medium);
            disp('assemble neighborhood coupling')
            obj = assemble_neighborhood_coupling_matrices(obj,ab5_table,positionArray,particleLmax,k_medium);            
        end
        
        function [y,obj] = couplingMatrixMultiply(obj,x)
            % return y = W*x
            
            y = x-x;
            
            % aggregation
            for jbox = 1:obj.numberOfBoxes
                xBox = x(obj.particleIdcs{jbox},:);
                xBox = xBox(:);
                obj.boxScatteredFieldCoefficients{jbox} = obj.particleAggregationMatrices{jbox} * xBox;
            end
            
            % box coupling
%             for jbox1 = 1:obj.numberOfBoxes %receiving box
%                 obj.boxIncomingFieldCoefficients{jbox1} = gpuArray(zeros(jmult_max(1,obj.lmax),1,'single'));
%             end
%             for idx = 1:length(obj.boxRelativeOffsetIdcs(:,1))
%                 idx
%                 mtr = gpuArray(squeeze(single(obj.boxCouplingMatrices(:,:,idx))));
%                 for jbox1 = 1:obj.numberOfBoxes %receiving box
%                     for jbox2 = 1:obj.numberOfBoxes %emitting box
%                         if all(obj.boxOffsetIdcs{jbox1} - obj.boxOffsetIdcs{jbox2} == obj.boxRelativeOffsetIdcs(idx,:))
%                             obj.boxIncomingFieldCoefficients{jbox1} = obj.boxIncomingFieldCoefficients{jbox1} + mtr*obj.boxScatteredFieldCoefficients{jbox2};
%                         end
%                     end
%                 end
%             end
                
                
            for jbox1 = 1:obj.numberOfBoxes %receiving box
                obj.boxIncomingFieldCoefficients{jbox1} = zeros(jmult_max(1,obj.lmax),1,'single');
                for jbox2 = 1:obj.numberOfBoxes %emitting box
                    if ~any(jbox2==obj.neighborBoxIdcs{jbox1})
                        relativeOffsetIdcs = obj.boxOffsetIdcs{jbox1} - obj.boxOffsetIdcs{jbox2};
                        [~,idx] = ismember(relativeOffsetIdcs,obj.boxRelativeOffsetIdcs,'rows');
                        obj.boxIncomingFieldCoefficients{jbox1} = obj.boxIncomingFieldCoefficients{jbox1} +obj.boxCouplingMatrices(:,:,idx)*obj.boxScatteredFieldCoefficients{jbox2};
                    end
                end
            end
            
            
            % disaggregation and neighborhood coupling
            for jbox = 1:obj.numberOfBoxes
                yBox = obj.particleDisaggregationMatrices{jbox} * obj.boxIncomingFieldCoefficients{jbox};
                for jbox2 = obj.neighborBoxIdcs{jbox}
                    xBoxNeighb = x(obj.particleIdcs{jbox2},:);
                    xBoxNeighb = xBoxNeighb(:);
                    ncm = gpuArray(obj.neighborhoodCouplingMatrices{jbox}{jbox2});
                    yBox = yBox + ncm * xBoxNeighb;
                end
                y(obj.particleIdcs{jbox},:) = reshape(yBox,length(obj.particleIdcs{jbox}),length(x(1,:)));
            end
        end
    end
end



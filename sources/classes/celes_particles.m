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

%> @file celes_particles.m
% ======================================================================
%> @brief Parameters that specify the particle aggregate
% ======================================================================

classdef celes_particles
    
    properties
        %> particle type, so far only 'sphere' implemented
        type  = 'sphere'  
        
        %> monodisperse or polydisperse? so far only 'mono' implemented. 
        %> that means that all particles are the same refractive index
        %> poly allows for differing refractive indices
        disperse  = 'mono'
        
        %> positions of the particles in the format [x(:),y(:),z(:)]
        positionArray
        
        %> complex refractive indices of the particles, n+ik
        refractiveIndexArray
        
        %> radii of the particles
        radiusArray
    end
    
    properties (Dependent)
        %> number of particles
        number
        
        %> maximal distance between two particles
        maxParticleDistance

        %> unique radii list
        uniqueRadii

        %> number of unique radii
        numUniqueRadii

        %> radiusArray in terms of indices given by uniqueRadii
        radiusArrayIndex
        
        %> unique index list
        uniqueRefractiveIndices
        
        %> number of unique refractive indices
        numUniqueRefractiveIndices
        
        %> refractiveIndexArray in terms of indices given by
        %  uniqueRefractiveIndices
        refractiveIndexArrayIndex
        
        %> unique pairs of refractive indices and radii
        %  for calculating mie coefficients
        uniqueRadiusIndexPairs
        
        %> unique pairs of refractive indices and radii
        %  for calculating mie coefficients
        uniqueSingleRadiusIndexPairs
        
        %> number of unique pairs of refractive indices and radii
        %  for calculating mie coefficients
        numUniquePairs
        
        %> single index encompassing radius and refractive index
        %  for indexing during matrix multiplication
        singleUniqueIndex
        
        %> singleUniqueIndexMap in terms of indices given by
        %  singleUniqueIndex
        %  serves as a lookup table for matrix multiplication
        singleUniqueArrayIndex
    end
    
    methods
        % ======================================================================
        %> @brief Set method for type
        % ======================================================================
        function obj = set.type(obj,value)
            switch value
                case 'sphere'
                    obj.type = value;
                otherwise
                    error('this particle type is at the moment not implemented')
            end
        end
        
        % ======================================================================
        %> @brief Set method for disperse
        % ======================================================================
        function obj = set.disperse(obj,value)
            switch value
                case 'mono'
                    obj.disperse=value;
                case 'poly'
                    obj.disperse=value;
                otherwise
                    error('this is at the moment not implemented')
            end
        end
        
        % ======================================================================
        %> @brief Set method for positionArray
        % ======================================================================
        function obj = set.positionArray(obj,value)
            if length(value(1,:))==3
                obj.positionArray = single(value);
            else
                error('illegal position array')
            end
        end
        
        % ======================================================================
        %> @brief Set method for refractive index
        % ======================================================================
        function obj = set.refractiveIndexArray(obj,value)
            obj.refractiveIndexArray=single(value);
        end
        
        % ======================================================================
        %> @brief Get method for unique refractive index values, returns ordered 
        %         vector of unique refractive indices
        % ======================================================================
        function value = get.uniqueRefractiveIndices(obj)
            value=unique(obj.refractiveIndexArray);
        end
        
        % ======================================================================
        %> @brief Get method for the number of unique refractive indices
        % ======================================================================
        function value = get.numUniqueRefractiveIndices(obj)
            value=length(unique(obj.refractiveIndexArray));
        end
        
        % ======================================================================
        %> @brief Get method for refractive index array in terms of indices given by 
        %        uniqueRefractiveIndices, sorted smallest to largest
        % ======================================================================
        function value = get.refractiveIndexArrayIndex(obj)
            value=dsearchn(obj.uniqueRefractiveIndices',obj.refractiveIndexArray');
        end        
        
        % ======================================================================
        %> @brief Set method for radiusArray
        %         added floor to guarantee validity of pairing function
        % ======================================================================
        function obj = set.radiusArray(obj,value)
            obj.radiusArray=single(value);
        end
        
        % ======================================================================
        %> @brief Get method for unique radii values, returns ordered vector of unique radii
        % ======================================================================
        function value = get.uniqueRadii(obj)
            value=unique(obj.radiusArray);
        end

        % ======================================================================
        %> @brief Get method for the number of unique radii
        % ======================================================================
        function value = get.numUniqueRadii(obj)
            value=length(unique(obj.radiusArray));
        end

        % ======================================================================
        %> @brief Get method for radius array in terms of indices given by uniqueRadii, sorted smallest to largest
        % ======================================================================
        function value = get.radiusArrayIndex(obj)
            value=dsearchn(obj.uniqueRadii',obj.radiusArray');
        end
        
        % ======================================================================
        %> @brief Get method to get unique pairs of radii and indices for
        %         computation of mie coefficients
        % ======================================================================
        function value = get.uniqueRadiusIndexPairs(obj)
            [radiiMap,indexMap] = meshgrid(obj.radiusArray,obj.refractiveIndexArray);
            allPairs = [radiiMap(:) indexMap(:)];
            value = unique(allPairs,'rows');
        end
        
        % ======================================================================
        %> @brief Get method to get unique pairs of radii and indices for
        %         matrix multiplication. For calculation of singleUniqueIndex
        % ======================================================================
        function value = get.uniqueSingleRadiusIndexPairs(obj)
            [radiiMap,indexMap] = meshgrid(obj.radiusArrayIndex,obj.refractiveIndexArrayIndex);
            allPairs = [radiiMap(:) indexMap(:)];
            value = unique(allPairs,'rows');
        end
        
        % ======================================================================
        %> @brief Get method for creating a single index encompassing unique pair 
        %         indices combined by pairing function:
        % p(rad,index) := 1/2(rad+index)(rad+index+1)+index
        % ======================================================================
        function value = get.singleUniqueIndex(obj)
            value = 1/2*(obj.uniqueSingleRadiusIndexPairs(:,1)+obj.uniqueSingleRadiusIndexPairs(:,2)).*(obj.uniqueSingleRadiusIndexPairs(:,1)+obj.uniqueSingleRadiusIndexPairs(:,2)+1)+obj.uniqueSingleRadiusIndexPairs(:,2);
        end
        
        % ======================================================================
        %> @brief Get method for a map of particles indexed by singleUniqueIndex 
        %         combined by pairing function:
        % p(rad,index) := 1/2(rad+index)(rad+index+1)+index
        % ======================================================================
        function value = get.singleUniqueArrayIndex(obj)
            pairedArray = 1/2*(obj.radiusArrayIndex+obj.refractiveIndexArrayIndex).*(obj.radiusArrayIndex+obj.refractiveIndexArrayIndex+1)+obj.refractiveIndexArrayIndex;
            value = dsearchn(obj.singleUniqueIndex,pairedArray);
        end
        
        % ======================================================================
        %> @brief Get method for the number of unique pairs of indices and
        %         radii
        % ======================================================================
        function value = get.numUniquePairs(obj)
            value = length(obj.uniqueRadiusIndexPairs(:,1));
        end
        
        % ======================================================================
        %> @brief Set method for particle number
        % ======================================================================
        function obj = set.number(obj,value)
            obj.positionArray=obj.positionArray(1:int32(value),:);
        end
        
        % ======================================================================
        %> @brief Set method for particle number
        % ======================================================================
        function value = get.number(obj)
            value=length(obj.positionArray(:,1));
        end
        
        % ======================================================================
        %> @brief Get method for particle number
        % ======================================================================
        function value = get.maxParticleDistance(obj)
            %value=max(pdist(obj.positionArray));  pdist part of statistics and machine learning toolbox and might not be available
            value=0;
            for jp1=1:obj.number
                diffs=bsxfun(@plus,obj.positionArray((jp1+1):end,:),-obj.positionArray(jp1,:));
                dists2 = diffs(:,1).^2+diffs(:,2).^2+diffs(:,3).^2;
                if max(dists2)>value^2
                    value=sqrt(max(dists2));
                end
            end
        end
        
        obj = sort_particles_by_partition(obj,partitioning);
        
        
    end
end


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
% ==============================================================================
%> @brief Parameters that specify the particle aggregate
% ==============================================================================

classdef celes_particles < matlab.System

    properties
        %> particle type, so far only 'sphere' implemented
        type = 'sphere'

        %> positions of the particles in the format [x(:),y(:),z(:)]
        positionArray

        %> complex refractive indices of the particles, n+ik
        refractiveIndexArray

        %> radii of the particles
        radiusArray
    end

    properties (SetAccess=private, Hidden)
        %> number of particles
        number

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
        %> uniqueRefractiveIndices
        refractiveIndexArrayIndex

        %> unique pairs of refractive indices and radii
        %> for calculating Mie coefficients
        uniqueRadiusIndexPairs

        %> unique pairs of refractive indices and radii
        %> for calculating Mie coefficients
        uniqueSingleRadiusIndexPairs

        %> number of unique pairs of refractive indices and radii
        %> for calculating Mie coefficients
        numUniquePairs

        %> single index encompassing radius and refractive index
        %> for indexing during matrix multiplication
        singleUniqueIndex

        %> singleUniqueIndexMap in terms of indices given by singleUniqueIndex
        %> serves as a lookup table for matrix multiplication
        singleUniqueArrayIndex

        %> maximal distance between two particles
        maxParticleDistance
    end

    methods
        % ======================================================================
        %> @brief Class constructor
        % ======================================================================
        function obj = celes_particles(varargin)
            setProperties(obj,nargin,varargin{:});
            obj.number = size(obj.positionArray,1);
            validatePropertiesImpl(obj);
            setupImpl(obj);
        end

        % ======================================================================
        %> @brief Compute unique values of refractive indices
        % ======================================================================
        function computeUniqueRefractiveIndices(obj)
            obj.uniqueRefractiveIndices = unique(obj.refractiveIndexArray);
            obj.numUniqueRefractiveIndices = length(obj.uniqueRefractiveIndices);
            obj.refractiveIndexArrayIndex = dsearchn(obj.uniqueRefractiveIndices, ...
                                                     obj.refractiveIndexArray);
        end

        % ======================================================================
        %> @brief Compute unique values of sphere radii
        % ======================================================================
        function computeUniqueRadii(obj)
            obj.uniqueRadii = unique(obj.radiusArray);
            obj.numUniqueRadii = length(obj.radiusArray);
            obj.radiusArrayIndex = dsearchn(obj.uniqueRadii, obj.radiusArray);
        end

        % ======================================================================
        %> @brief Compute unique combinations of radii and refractive indices
        % ======================================================================
        function computeUniqueRadiiIndexPairs(obj)
            obj.uniqueRadiusIndexPairs = ...
                unique([obj.radiusArray, obj.refractiveIndexArray],'rows');
            obj.uniqueSingleRadiusIndexPairs = ...
                unique([obj.radiusArrayIndex, obj.refractiveIndexArrayIndex],'rows');
            obj.radiusArrayIndex = dsearchn(obj.uniqueRadii, obj.radiusArray);
        end

        % ======================================================================
        %> @brief Compute single unique index
        % ======================================================================
        function computeSingleUniqueIndex(obj)
            obj.singleUniqueIndex = 1/2*(obj.uniqueSingleRadiusIndexPairs(:,1)+ ...
                                         obj.uniqueSingleRadiusIndexPairs(:,2)).* ...
                                        (obj.uniqueSingleRadiusIndexPairs(:,1)+ ...
                                         obj.uniqueSingleRadiusIndexPairs(:,2)+1)+ ...
                                         obj.uniqueSingleRadiusIndexPairs(:,2);
            pairedArray = 1/2*(obj.radiusArrayIndex+ ...
                               obj.refractiveIndexArrayIndex).* ...
                              (obj.radiusArrayIndex+ ...
                               obj.refractiveIndexArrayIndex+1)+ ...
                               obj.refractiveIndexArrayIndex;
            obj.singleUniqueArrayIndex = dsearchn(obj.singleUniqueIndex, pairedArray);
            obj.numUniquePairs = length(obj.uniqueRadiusIndexPairs(:,1));
        end

        % ======================================================================
        %> @brief Set the maximalParticleDistance attribute to the correct value
        % ======================================================================
        function compute_maximal_particle_distance(obj)
            try
                pos = convhull(double(obj.positionArray)); % double required by convhull
                pos = obj.positionArray(unique(pos(:)),:); % unique vertices
            catch % convhull fails if particles are coplanar or collinear
                pos = obj.positionArray;
            end
            try
                obj.maxParticleDistance = max(pdist(pos));
            catch % in case pdist is not available
                % drop-in replacement: a bit slower, larger memory footprint
                % obj.maxParticleDistance = sqrt(max(max(sum(bsxfun(@minus,reshape(size(pos,[size(pos,1),1,3]),reshape(pos,[1,size(pos,1),3])).^2,3))));
                %
                % quick&dirty overestimate, probably OK in most (all?) cases
                obj.maxParticleDistance = sqrt(sum((min(pos)-max(pos)).^2));
            end
        end
    end

    methods(Access = protected)
        % ======================================================================
        %> @brief Class implementation
        % ======================================================================
        function setupImpl(obj)
            computeUniqueRefractiveIndices(obj)
            computeUniqueRadii(obj)
            computeUniqueRadiiIndexPairs(obj)
            computeSingleUniqueIndex(obj)
            compute_maximal_particle_distance(obj);
        end

        % ======================================================================
        %> @brief Validate class properties
        % ======================================================================
        function validatePropertiesImpl(obj)
            if lower(obj.type) ~= 'sphere'
                error('this particle type is not implemented at the moment')
            end
            try validateattributes(obj.positionArray,{'numeric'},{'real','nonnan','finite','2d','ncols',3})
            catch e, error('invalid positionArray: %s', e.message); end
            try validateattributes(obj.refractiveIndexArray,{'numeric'},{'nonnan','finite','numel',obj.number})
            catch e, error('invalid refractiveIndexArray: %s', e.message); end
            try validateattributes(obj.radiusArray,{'numeric'},{'nonnan','real','finite','numel',obj.number})
            catch e, error('invalid radiusArray: %s', e.message); end
        end
    end
end

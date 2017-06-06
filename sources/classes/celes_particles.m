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
        %> that means that all particles are the same.
        disperse  = 'mono'
        
        %> positions of the particles in the format [x(:),y(:),z(:)]
        positionArray
        
        %> complex refractive indices of the particles, n+ik
        refractiveIndex
        
        %> radius of the particles
        radius
        
        %> maximal distance between two particles
        maxParticleDistance
    end
    
    properties (Dependent)
        %> number of particles
        number
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
        function obj = set.refractiveIndex(obj,value)
            obj.refractiveIndex=single(value);
        end
        
        % ======================================================================
        %> @brief Set method for radius
        % ======================================================================
        function obj=set.radius(obj,value)
            obj.radius=single(value);
        end
        
        % ======================================================================
        %> @brief Set method for particle number
        % ======================================================================
        function obj=set.number(obj,value)
            obj.positionArray=obj.positionArray(1:int32(value),:);
        end
        
        % ======================================================================
        %> @brief Set method for particle number
        % ======================================================================
        function value=get.number(obj)
            value=length(obj.positionArray(:,1));
        end
        
        % ======================================================================
        %> @brief Set the maximalParticleDistance attribute to the correct value
        % ======================================================================
        function obj = compute_maximal_particle_distance(obj)
            %value=max(pdist(obj.positionArray));  pdist part of statistics and machine learning toolbox and might not be available
            obj.maxParticleDistance=0;
            for jp1=1:obj.number
                diffs=bsxfun(@plus,obj.positionArray((jp1+1):end,:),-obj.positionArray(jp1,:));
                dists2 = diffs(:,1).^2+diffs(:,2).^2+diffs(:,3).^2;
                if max(dists2)>obj.maxParticleDistance^2
                    obj.maxParticleDistance=sqrt(max(dists2));
                end
            end
        end
        
        obj = sort_particles_by_partition(obj,partitioning);
        
        
    end
end


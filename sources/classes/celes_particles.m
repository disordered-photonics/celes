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
    end
    
    properties (Dependent)
        %> number of particles
        number
        
        %> maximal distance between two particles
        maxParticleDistance
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
        %> @brief Get method for particle number
        % ======================================================================
        function value=get.maxParticleDistance(obj)
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


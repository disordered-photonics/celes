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

%> @file celes_simulation.m
% ======================================================================
%> @brief Central data structure of the celes software
%
%> The simulation class contains all input, intermediate results and 
%> output for one calculation.
% ======================================================================
classdef celes_simulation
    
    properties
        %> celes_input object which contains the parameters that specify
        %> the simulation geometry and initial field
        input
        
        %> celes_numerics object which contains the numerical settings
        numerics
        
        %> celes_tables object which contains lookup tables and other
        %> intermediate results
        tables
        
        %> celes_output object which contains the results of the
        %> simulation
        output

    end
    
    properties (Dependent)
        %> single array which contains a grid of distances used for the
        %> lookup of the spherical hankel function in the particle coupling
        lookupParticleDistances
       
    end
    
    methods
        % ======================================================================
        %> @brief get method for dependent property lookupParticleDistances
        % ======================================================================
        function value = get.lookupParticleDistances(obj)
            value = [0,0:obj.numerics.particleDistanceResolution:(obj.input.particles.maxParticleDistance + 3*obj.numerics.particleDistanceResolution)];  % add two zeros at beginning to allow for interpolation also in first segment
        end

        % ======================================================================
        %> @brief Evaluate the Mie coefficients
        %>
        %> @return celes_simulation object with updated mieCoefficients
        % ======================================================================
        function obj = computeMieCoefficients(obj)
            fprintf(1,'compute Mie coefficients ...');
            switch obj.input.particles.type
                case 'sphere'
                    obj.tables.mieCoefficients = zeros(obj.input.particles.numUniqueRadii,obj.numerics.nmax,'single');
                for r_i=1:obj.input.particles.numUniqueRadii
                    for tau=1:2
                         for l=1:obj.numerics.lmax
                             for m=-l:l
                                jmult = multi2single_index(1,tau,l,m,obj.numerics.lmax);
                   	            obj.tables.mieCoefficients(r_i,jmult) = T_entry(tau,l,obj.input.k_medium,obj.input.k_particle,obj.input.particles.uniqueRadii(r_i));
                             end
                         end
                    end
                end
                otherwise
                    error('particle type not implemented')
            end
            fprintf(1,' done\n');
        end
        
        % ======================================================================
        %> @brief Prepare a lookup for the a and b coefficients for
        %> particle coupling
        %>
        %> @return celes_simulation object with updated translationTable
        % ======================================================================
        function obj = computeTranslationTable(obj)
            fprintf(1,'compute translation table ...');
            obj.tables.translationTable = translation_table_ab(obj.numerics.lmax);
            fprintf(1,' done\n');
        end
        
        % ======================================================================
        %> @brief Evaluate the initial field coefficients \f$a^S_{0,n}\f$
        %> of the initial field expansion around each particle:
        %> \f$\mathbf{E}_0=\sum_n a^S_{0,n}\mathbf{\Psi}^{(1)}_n\f$
        %>
        %> @return celes_simulation object with updated initialFieldCoefficients
        % ======================================================================
        function obj = computeInitialFieldCoefficients(obj)
            fprintf(1,'compute initial field coefficients ...');
            if isfinite(obj.input.initialField.beamWidth) && obj.input.initialField.beamWidth
                fprintf(1,' Gaussian beam ...');
                if obj.input.initialField.normalIncidence
                    obj.tables.initialFieldCoefficients = initial_field_coefficients_wavebundle_normal_incidence(obj);
                else
                    error('this case is not implemented')
                end
            else % infinite or 0 beam width
                fprintf(1,' plane wave ...');
                obj.tables.initialFieldCoefficients = initial_field_coefficients_planewave(obj);
            end
            fprintf(1,' done\n');
        end
        
        % ======================================================================
        %> @brief Evaluate the power flux of the initial field
        %>
        %> @return celes_simulation object with updated initialFieldPower
        % ======================================================================
        function obj = computeInitialFieldPower(obj)
            fprintf(1,'compute initial field power ...');
            if obj.input.initialField.normalIncidence
                obj.output.initialFieldPower = initial_power_wavebundle_normal_incidence(obj);
            else
                error('this case is not implemented')
            end
            fprintf(1,' done\n');
        end
        
        % ======================================================================
        %> @brief Compute the scattered field coefficients b by iteratively 
        %> solving the linear system M*b=T*aI
        %>
        %> @param Optional: b0, initial guess for scattered field
        %> coefficients
        %> @return celes_simulation object with updated initialFieldPower
        % ======================================================================
        function obj = computeScatteredFieldCoefficients(obj,varargin)
            fprintf(1,'compute scattered field coefficients ...');
            mmm = @(x) obj.masterMatrixMultiply(x);
            [b,convHist] = obj.numerics.solver.run(mmm,obj.tables.rightHandSide(:),varargin{:});
            obj.tables.scatteredFieldCoefficients = reshape(gather(b),size(obj.tables.rightHandSide));
            obj.output.convergenceHistory = convHist;
        end
        
        % ======================================================================
        %> @brief Compute the plane wave pattern of the scattered field
        %> (i.e., the expansion coefficients of the scattered field in plane
        %> vector wave functions)
        %>
        %> @return celes_simulation object with updated
        %> output.scatteredFieldPlaneWavePattern
        % ======================================================================
        function obj = computeScatteredFieldPWP(obj)
            fprintf(1,'compute scattered field plane wave pattern: ');
            obj.output.scatteredFieldPlaneWavePattern = scattered_field_plane_wave_pattern(obj);
            fprintf(1,' ... done\n');
        end
        
        % ======================================================================
        %> @brief Compute the plane wave pattern of the total field
        %> (i.e., the expansion coefficients of the scattered field in plane
        %> vector wave functions)
        %>
        %> @return celes_simulation object with updated
        %> output.totalFieldPlaneWavePattern
        % ======================================================================
        function obj = computeTotalFieldPWP(obj)
            fprintf(1,'compute total field coefficients table ...');
            obj.output.totalFieldPlaneWavePattern=cell(2,1);
            pwpScat = obj.output.scatteredFieldPlaneWavePattern;
            pwpIn = initial_field_plane_wave_pattern(obj);
            for pol=1:2
                obj.output.totalFieldPlaneWavePattern{pol} = pwpScat{pol};
                obj.output.totalFieldPlaneWavePattern{pol} = obj.output.totalFieldPlaneWavePattern{pol}.addTo(pwpIn{pol});
            end
            fprintf(1,' done\n');
        end
        
        % ======================================================================
        %> @brief Evaluate the power flux of the total field, both in
        %> forward and in backward direction
        %>
        %> @return celes_simulation object with updated
        %> output.totalFieldForwardPower and output.totalFieldBackwardPower
        % ======================================================================
        function obj = computeTotalFieldPower(obj)
            fprintf(1,'compute total field power ...');
            obj.output.totalFieldForwardPower = gather(pwp_power_flux(obj.output.totalFieldPlaneWavePattern{1},obj,'forward') + pwp_power_flux(obj.output.totalFieldPlaneWavePattern{2},obj,'forward'));
            obj.output.totalFieldBackwardPower = gather(pwp_power_flux(obj.output.totalFieldPlaneWavePattern{1},obj,'backward') + pwp_power_flux(obj.output.totalFieldPlaneWavePattern{2},obj,'backward'));
            fprintf(1,' done\n');
        end
        
        % ======================================================================
        %> @brief First prepare the scattered and total field's plane wave
        %> pattern, then evaluate the power flux
        %>
        %> @return celes_simulation object with updated
        %> scatteredFieldPlaneWavePattern, totalFieldPlaneWavePattern,
        %> output.totalFieldForwardPower and output.totalFieldBackwardPower
        % ======================================================================
        function obj = evaluatePower(obj)
            tpow = tic;
            obj = obj.computeScatteredFieldPWP;
            obj = obj.computeTotalFieldPWP;
            obj = obj.computeTotalFieldPower;
            obj.output.powerEvaluationTime = toc(tpow);
            fprintf(1,'power flux evaluated in %.1f seconds.\n',obj.output.powerEvaluationTime);
        end
        
        % ======================================================================
        %> @brief Evaluate the initial (near)field at the positions 
        %> specified in the input. The field can then be plotted.
        %>
        %> @return celes_simulation object with updated 
        %> output.InitialField
        % ======================================================================
        function obj = evaluateInitialField(obj)
            fprintf(1,'evaluate initial field ...');
            obj.output.initialField = compute_initial_field(obj);
            fprintf(1,' done\n');
        end
        
        % ======================================================================
        %> @brief Evaluate the scattered (near)field at the positions 
        %> specified in the input. The field can then be plotted.
        %>
        %> @return celes_simulation object with updated 
        %> output.scatteredField
        % ======================================================================
        function obj = evaluateScatteredField(obj)
            fprintf(1,'evaluate scattered field ...');
            obj.output.scatteredField = compute_scattered_field(obj);
            fprintf(1,' done\n');
        end
        
        % ======================================================================
        %> @brief Evaluate the internal (near)field at the positions 
        %> specified in the input. The field can then be plotted.
        %>
        %> @return celes_simulation object with updated 
        %> output.internalField
        % ======================================================================
        function obj = evaluateInternalField(obj)
            fprintf(1,'evaluate internal field ...');
            [obj.output.internalField,obj.output.internalIndices] = compute_internal_field(obj);
            fprintf(1,' done\n');
        end
        
        % ======================================================================
        %> @brief Evaluate both the initial and the scattered (near)field 
        %> at the positions specified in the input. The field can then be
        %> plotted.
        %>
        %> @return celes_simulation object with updated 
        %> output.initialField and output.scatteredField
        % ======================================================================
        function obj = evaluateFields(obj)
            tfld = tic;
            obj = obj.evaluateInitialField;
            obj = obj.evaluateScatteredField;
            obj = obj.evaluateInternalField;
            obj.output.fieldEvaluationTime = toc(tfld);
            fprintf(1,'fields evaluated in %.1f seconds.\n',obj.output.fieldEvaluationTime);
        end
        
        % ======================================================================
        %> @brief Multiply the master matrix M=1-T*W to some vector x
        %>
        %> @param Vector x of incoming field SVWF coefficients
        %> @return Vector M*x
        % ======================================================================
        function Mx = masterMatrixMultiply(obj,value)
            value=value(:);
            Wx=coupling_matrix_multiply(obj,value);
            Wx=reshape(Wx,obj.input.particles.number,obj.numerics.nmax);
            TWx = obj.tables.mieCoefficients(obj.input.particles.radiusArrayIndex,:).*Wx;
            Mx = value - TWx(:);
        end
        
        % ======================================================================
        %> @brief Run the simulation.
        %> 
        %> A simulation run includes:
        %> - computation of initial field power
        %> - computation of Mie coefficients
        %> - computation of the translation table
        %> - preparation of the particle partitioning (if blockdiagonal
        %>   preconditioner is active)
        %> - preparation of the blockdiagonal preconditioner (if active)
        %> - computation of initial field coefficients
        %> - solution of linear system
        %> 
        %> @param Optional: Initial guess b0 for the scattered field
        %> coefficients vector b
        %> @return celes_simulation object with various fields updated
        % ======================================================================
        function obj = run(obj,varargin)
            print_logo
            print_parameters(obj)
            tcomp=tic;
            cuda_compile(obj.numerics.lmax);
            fprintf(1,'starting simulation.\n');
            obj = obj.computeInitialFieldPower;
            obj = obj.computeMieCoefficients; %modify for radius
            obj = obj.computeTranslationTable;
            tprec=tic;
            if strcmp(obj.numerics.solver.preconditioner.type,'blockdiagonal')
                fprintf(1,'make particle partition ...');
                partitioning = make_particle_partion(obj.input.particles.positionArray,obj.numerics.solver.preconditioner.partitionEdgeSizes);
                if obj.input.particles.numUniqueRadii == 1
                    obj = sort_particles_by_partition(obj,partitioning);
                    partitioning = make_particle_partion(obj.input.particles.positionArray,obj.numerics.solver.preconditioner.partitionEdgeSizes);
                end
                obj.numerics.solver.preconditioner.partitioning = partitioning;
                fprintf(1,' done\n');
                obj = obj.numerics.solver.preconditioner.prepare(obj);
            end
            obj.output.preconiditionerPreparationTime = toc(tprec);
            fprintf(1,'preconditioner prepared in %.1f seconds.\n',obj.output.preconiditionerPreparationTime);
            tsolv=tic;
            obj = obj.computeInitialFieldCoefficients;
            obj = obj.computeScatteredFieldCoefficients(varargin{:});
            obj.output.solverTime = toc(tsolv);
            fprintf(1,'solver terminated in %.1f seconds.\n',obj.output.solverTime);
            obj.numerics.solver.preconditioner.factorizedMasterMatrices = []; % clear memory intensive fields
            obj.numerics.solver.preconditioner.masterMatrices = [];
            obj.output.runningTime = toc(tcomp);
            fprintf(1,'simulation ran in %.1f seconds.\n',obj.output.runningTime);
        end
    end
end

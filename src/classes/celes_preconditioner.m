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

%> @file celes_preconditioner.m
% ==============================================================================
%> @brief Preconditioner for the iterative celes_solver objects
% ==============================================================================

classdef celes_preconditioner < matlab.System

    properties
        %> what kind of preconditioner?
        %> at the moment, only 'none' and 'blockdiagonal' implemented
        type = 'none'

        %> for type blockdiagonal: edge sizes of the cuboids that define
        %> the partition of the particles (i.e. the blocks)
        %> format: [dx,dy,dz]
        partitionEdgeSizes

        %> for type blockdiagonal: cell array of particle indices that are
        %> grouped into the same partition
        partitioning

        %> for type blockdiagonal: cell array of SVWF coefficient indices
        %> that belong to the corresponding partitioning
        partitioningIdcs

        %> for type blockdiagonal: cell array of blocks of master matrix
        %> M=1-TW that correspond to one partitioning
        masterMatrices

        %> for type blockdiagonal: cell array of LU decompositions of
        %> masterMatrices
        factorizedMasterMatrices
    end

    methods
        % ======================================================================
        %> @brief Class constructor
        % ======================================================================
        function obj = celes_preconditioner(varargin)
            setProperties(obj,nargin,varargin{:});
            validatePropertiesImpl(obj);
        end

        % ======================================================================
        %> @brief Run the preconditioner.
        %>
        %> @param Coefficient vector x
        %> @return y = M'\x
        %> where M' is a blockdiagonal approximation to the linear system M
        % ======================================================================
        function simul = prepare(obj,simul)
            switch lower(obj.type)
                case 'blockdiagonal'
                    fprintf(1,'prepare blockdiagonal preconditioner ...\n');
                    msg = '';
                    lmax = simul.numerics.lmax;
                    nmax = simul.numerics.nmax;
                    k = simul.input.k_medium;

                    for jp = 1:length(obj.partitioning)
                        spherArr = obj.partitioning{jp};
                        NSi = length(spherArr);

                        Idcs= [];
                        for n=1:nmax
                            Idcs = [Idcs; ...
                                    spherArr+ ...
                                    simul.input.particles.number*(n-1)];
                        end
                        simul.numerics.solver.preconditioner.partitioningIdcs{jp} = Idcs;

                        fprintf(1,repmat('\b',[1,length(msg)]));
                        msg = sprintf('partition %i of %i: %i particles -- compute master matrix ...',jp,length(obj.partitioning),length(spherArr));
                        fprintf(1,msg);

                        M = eye(NSi*simul.numerics.nmax,'single');

                        [x2,x1] = meshgrid(simul.input.particles.positionArray(spherArr,1));
                        [y2,y1] = meshgrid(simul.input.particles.positionArray(spherArr,2));
                        [z2,z1] = meshgrid(simul.input.particles.positionArray(spherArr,3));
                        x1mnx2 = x1-x2;
                        y1mny2 = y1-y2;
                        z1mnz2 = z1-z2;
                        dTab = sqrt(x1mnx2.^2+y1mny2.^2+z1mnz2.^2);
                        ctTab = z1mnz2./dTab;
                        stTab = sqrt(1-ctTab.^2);
                        phiTab = atan2(y1mny2,x1mnx2);
                        Plm = legendre_normalized_trigon(ctTab,stTab,2*lmax);
                        particleArrayInd = simul.input.particles.singleUniqueArrayIndex;
                        for p = 0:2*lmax
                            sphHank = sph_bessel(3,p,k*dTab);
                            for tau1 = 1:2
                                for l1 = 1:lmax
                                    for m1 = -l1:l1
                                        n1 = multi2single_index(1,tau1,l1,m1,lmax);
                                        n1S1Arr = (1:NSi)+(n1-1)*NSi;
                                        for tau2 = 1:2
                                            for l2 = 1:lmax
                                                for m2 = -l2:l2
                                                    if abs(m1-m2) <= p
                                                        n2 = multi2single_index(1,tau2,l2,m2,lmax);
                                                        n2S2Arr = (1:NSi)+(n2-1)*NSi;
                                                        TWn1n2 = simul.tables.mieCoefficients(particleArrayInd(spherArr(:)),n1)*simul.tables.translationTable.ab5(n2,n1,p+1) .* Plm{p+1,abs(m1-m2)+1} .* sphHank .* exp(1i*(m2-m1)*phiTab) ;
                                                        s1eqs2 = logical(eye(NSi));
                                                        TWn1n2(s1eqs2(:)) = 0; % jS1=jS2
                                                        M(n1S1Arr,n2S2Arr) = M(n1S1Arr,n2S2Arr) - TWn1n2;
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end

                        fprintf(1,repmat('\b',[1,length(msg)]));
                        msg = sprintf('partition %i of %i: %i particles -- factorize master matrix ...',jp,length(obj.partitioning),length(spherArr));
                        fprintf(1,msg);

                        [Y,U,P] = lu(simul.numerics.deviceArray(M),'vector');
                        Y(U~=0) = U(U~=0);
                        simul.numerics.solver.preconditioner.masterMatrices{jp} = M;
                        simul.numerics.solver.preconditioner.factorizedMasterMatrices{jp}.Y = gather(Y);
                        simul.numerics.solver.preconditioner.factorizedMasterMatrices{jp}.P = gather(P);
                    end
                    fprintf(' done\n')
                case 'none'
            end
        end

        % ======================================================================
        %> @brief Run the preconditioner.
        %>
        %> @param Coefficient vector x
        %> @return y = M'\x
        %> where M' is a blockdiagonal approximation to the linear system M
        % ======================================================================
        function value = run(obj,rhs,varargin)

            switch lower(obj.type)
                case 'blockdiagonal'

                    if isempty(varargin)
                        verbose = false;
                    else
                        verbose = varargin{1};
                    end
                    if verbose
                        fprintf('apply preconditioner ... ')
                    end
                    prec_timer = tic;

                    optL.LT = true; % settings for linsolve
                    optU.UT = true;
                    rhs = gather(rhs(:));
                    value = zeros(size(rhs),'like',rhs);
                    for jp = 1:length(obj.partitioning)
                        p = obj.factorizedMasterMatrices{jp}.P;
                        U = triu(obj.factorizedMasterMatrices{jp}.Y);
                        L = tril(obj.factorizedMasterMatrices{jp}.Y,-1)+eye(length(p),'single');
                        value(obj.partitioningIdcs{jp}) = linsolve (U, linsolve (L, rhs(obj.partitioningIdcs{jp}(p)), optL), optU);
                    end
                    value = gpuArray(value(:));
                    prec_time = toc(prec_timer);
                    if verbose
                        fprintf('done in %f seconds.\n', prec_time)
                    end

                case 'none'
                    value=rhs;
            end
        end
    end

    methods(Access = protected)
        % ======================================================================
        %> @brief Validate class properties
        % ======================================================================
        function validatePropertiesImpl(obj)
            try validateattributes(obj.type,{'char'},{'nonempty'})
            catch e, error('invalid preconditioner type: %s', e.message); end
            if strcmpi(obj.type,'blockdiagonal')
                try validateattributes(obj.partitionEdgeSizes,{'numeric'},{'real','nonnan','finite','row','numel',3})
                catch e, error('invalid partitionEdgeSizes array: %s', e.message); end
            end
        end
    end
end

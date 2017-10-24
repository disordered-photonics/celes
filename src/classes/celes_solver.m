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

%> @file celes_solver.m
% ==============================================================================
%> @brief Solver settings for the master equation
%>        \f$\sum_{i',n'}M^{ii'}_{nn'} b_{n'}^{i'} = \sum_{n'} T^i_{nn'} a_{\mathrm{in},n'}^i\f$
% ==============================================================================

classdef celes_solver < matlab.System

    properties
        %> solver type, at the moment 'BiCGStab' and 'GMRES' implemented
        type = 'gmres';

        %> solution tolerance. solver stops when relative error smaller
        tolerance = 1e-4;

        %> maximal number of iterations. solver stops if exceeded
        maxIter = 1000;

        %> for type='GMRES': restart parameter
        restart = 100;

        %> celes_preconditioner object for faster convergence of solver
        preconditioner = celes_preconditioner;
    end

    methods
        % ======================================================================
        %> @brief Class constructor
        % ======================================================================
        function obj = celes_solver(varargin)
            setProperties(obj,nargin,varargin{:});
            validatePropertiesImpl(obj);
        end

        % ======================================================================
        %> @brief Run the solver
        %>
        %> @param mmm (function handle): master matrix multiply operator M*
        %> @param rhs (float array): right hand side of linear system
        %> @param Optional: y0 (float array): initial guess for solution vector
        %> @param Optional: verbose (logical): If true (default), display
        %>        detailed timing information
        %> @return y = M\rhs, an approximation to the solution of the linear system
        % ======================================================================
        function [value,convergenceHistory] = run(obj,mmm,rhs,varargin)
            if isempty(varargin)
                initial_guess = gather(rhs);
            else
                initial_guess = varargin{1};
            end
            if length(varargin) > 1
                verbose = varargin{2};
            else
                verbose = false;
            end

            prh = @(x) gather(obj.preconditioner.run(x,verbose));

            switch lower(obj.type)
                case 'bicgstab'
                    [value,~,~,~,convergenceHistory] = ...
                        bicgstab_custom(mmm,rhs,obj.tolerance,obj.maxIter,prh,[],initial_guess);
                case 'gmres'
                    fh = str2func('gmres');
                    try
                        fetch_and_patch_gmres();
                        fh = str2func('gmres_monitor');
                    catch
                        warning('using MATLAB default GMRES');
                    end
                    [value,~,~,~,convergenceHistory] = ...
                        fh(mmm,rhs,obj.restart,obj.tolerance,obj.maxIter,prh,[],initial_guess);
                otherwise
                    error('please specify a valid solver type')
            end
        end
    end

    methods(Access = protected)
        % ======================================================================
        %> @brief Validate class properties
        % ======================================================================
        function validatePropertiesImpl(obj)
            try validateattributes(obj.type,{'char'},{'nonempty'})
            catch e, error('invalid solver type: %s', e.message); end
            try validateattributes(obj.tolerance,{'numeric'},{'real','nonnan','finite','nonempty','positive','scalar'})
            catch e, error('invalid tolerance value: %s', e.message); end
            try validateattributes(obj.maxIter,{'numeric'},{'nonempty','integer','positive'})
            catch e, error('invalid maxIter value: %s', e.message); end
            try validateattributes(obj.restart,{'numeric'},{'nonempty','integer','positive'})
            catch e, error('invalid restart value: %s', e.message); end
        end
    end
end

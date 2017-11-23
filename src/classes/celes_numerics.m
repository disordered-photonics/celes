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

%> @file celes_numerics.m
% ==============================================================================
%> @brief Parameters describing the numerical settings used in the simulation
% ==============================================================================

classdef celes_numerics < matlab.System

    properties
        %> maximal polar multipole order
        lmax

        %> array of polar angles (between 0 and pi) that specifies the
        %> discretization of angles in the plane wave patterns
        polarAnglesArray

        %> array of azimuthal angles (between 0 and 2*pi) that specifies
        %> the discretization of angles in the plane wave patterns
        azimuthalAnglesArray

        %> shall the GPU be used in computations of, e.g., the
        %> preconditioner?
        %> NOTE: In the current version, the GPU is always used for the
        %> matrix vector products. Future versions might offer the
        %> opportunity to run the whole calculation on the CPU
        gpuFlag = true

        %> celes_solver object that contains all the information about
        %> how the linear system is to be solved
        solver = celes_solver

        %> resolution of the radial particle distance in the lookup tables
        %> of the translation operator
        particleDistanceResolution = single(10);
    end

    properties (SetAccess=private, Hidden)
        %> number of unknowns per particle
        nmax

        %> single precision (lmax+1)x(lmax+1)x(ceil(lmax/2)) array of
        %> coefficients of the (trigonometric) associated legendre polynomials
        Plm_coeff_table
    end

    methods
        % ======================================================================
        %> @brief Class constructor
        % ======================================================================
        function obj = celes_numerics(varargin)
            setProperties(obj,nargin,varargin{:});
            validatePropertiesImpl(obj);
            setupImpl(obj);
        end

        % ======================================================================
        %> @brief Compute nmax
        % ======================================================================
        function computeNmax(obj)
            obj.nmax = jmult_max(1,obj.lmax);
        end

        % ======================================================================
        %> @brief Compute the coefficients for the calculation of the associated
        %>        Legendre functions (requires MATLAB's symbolic toolbox)
        %>
        %> @param lmax
        %> @retval Plm_coeff_table (single precision float array of size
        %>         lmax+1 x lmax+1 x ceil(lmax/2)). Plm_coeff_table can be used
        %>         to reconstruct the actual Legendre functions for example by
        %>         using the function check_legendre()
        % ======================================================================
        function Plm_coefficients(obj)
            syms ct
            syms st
            plm = legendre_normalized_trigon(ct,st,2*obj.lmax);

            obj.Plm_coeff_table = zeros(2*obj.lmax+1, ...
                                        2*obj.lmax+1, ...
                                        ceil(2*obj.lmax/2), ...
                                        'single');

            for l = 0:2*obj.lmax
                for m = 0:l
                    [cf,~] = coeffs(plm{l+1,m+1});
                    obj.Plm_coeff_table(l+1,m+1,1:length(cf)) = single(cf);
                end
            end
%             clear ct st
%             reset(symengine)
        end

        % ======================================================================
        %> @brief Returns a gpuArray if gpuFlag, and a usual array
        %> otherwise
        %>
        %> @param Array x
        %> @return gpuArray(x) if gpuFlag, otherwise x
        % ======================================================================
        function arr_out = deviceArray(obj,arr_in)
            if obj.gpuFlag
                arr_out = gpuArray(arr_in);
            else
                arr_out = arr_in;
            end
        end

    end

    methods(Access = protected)
        % ======================================================================
        %> @brief Class implementation
        % ======================================================================
        function setupImpl(obj)
            computeNmax(obj);
            Plm_coefficients(obj);
        end

        % ======================================================================
        %> @brief Validate class properties
        % ======================================================================
        function validatePropertiesImpl(obj)
            try validateattributes(obj.lmax,{'numeric'},{'positive','integer','scalar'})
            catch e, error('invalid lmax value: %s', e.message); end
            try validateattributes(obj.polarAnglesArray,{'numeric'},{'real','nonnan','increasing','row','>=',0,'<=',pi})
            catch e, error('invalid polarAnglesArray array: %s', e.message); end
            try validateattributes(obj.azimuthalAnglesArray,{'numeric'},{'real','nonnan','increasing','row','>=',0,'<=',2*pi})
            catch e, error('invalid azimuthalAnglesArray array: %s', e.message); end
            try validateattributes(obj.particleDistanceResolution,{'numeric'},{'real','nonnan','positive','scalar'})
            catch e, error('invalid particleDistanceResolution value: %s', e.message); end
            try validateattributes(obj.solver,{'celes_solver'},{'nonempty'})
            catch, error('expected a valid celes_solver instance'); end
        end
    end
end

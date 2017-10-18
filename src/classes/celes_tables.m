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

%> @file celes_tables.m
% ==============================================================================
%> @brief Objects of this class hold large tables of interim results
% ==============================================================================

classdef celes_tables

    properties
        %> a table of coefficients needed for the SVWF translation
        translationTable

        %> T-matrix of the spheres
        mieCoefficients

        %> index array to map T-matrices to particles
        singleUniqueArrayIndex

        %> coefficients of the regular SVWF expansion of the initial
        %> excitation
        initialFieldCoefficients

        %> coefficients of the outgoing SVWF expansion of the scattered
        %> field
        scatteredFieldCoefficients
    end

    properties (Dependent)
        %> right hand side T*aI of linear system M*b=T*aI
        rightHandSide
    end

    methods
        % ======================================================================
        %> @brief Get method for rightHandSide
        % ======================================================================
        function TaI = get.rightHandSide(obj)
            TaI = obj.mieCoefficients(obj.singleUniqueArrayIndex,:).* ...
                  obj.initialFieldCoefficients;
        end
    end
end

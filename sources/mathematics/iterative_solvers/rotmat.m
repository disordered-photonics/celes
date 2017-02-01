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

%======================================================================
%> @brief Givens rotation matrix - helper function for gmres_custom()
%======================================================================
function [ c, s ] = rotmat( a, b )

% from the package "Templates for the Solution of Linear Systems"
% http://it.mathworks.com/matlabcentral/fileexchange/2158-templates-for-the-solution-of-linear-systems

% The authors of the original routine provided the following copyright
% statement:
%
% "All M-files are copyrighted, 1993, by Richard Barrett, Michael Berry,
% Tony Chan, James Demmel, June Donato, Jack Dongarra, Victor Eijkhout,
% Roldan Pozo, Charles Romine, and Henk van der Vorst.
%
% These M-files are User Contributed Routines that have been contributed to
% The MathWorks, Inc. and which are being redistributed by The MathWorks, upon
% request, on an "as is" basis.  A User Contributed Routine is not a product of
% The MathWorks and The MathWorks assumes no responsibility for any errors
% that may exist in the routines.
%
% The M-files were created to supplement "Templates for the Solution of
% Linear Systems:  Building Blocks for Iterative Methods," by Richard Barrett,
% Michael Berry, Tony Chan, James Demmel, June Donato, Jack Dongarra,
% Victor Eijkhout, Roldan Pozo, Charles Romine, and Henk van der Vorst (SIAM, 1994).
% You are free to modify any of the files and create new functions, provided
% that you acknowledge the source in any publication and do not sell the modified file.
%
% The M-files have been created using MATLAB Version 4.0 and have been
% written and checked for use on The Student Edition of MATLAB, 1992."

%
% Compute the Givens rotation matrix parameters for a and b.
%
   if ( b == 0.0 ),
      c = 1.0;
      s = 0.0;
   elseif ( abs(b) > abs(a) ),
      temp = a / b;
      s = 1.0 / sqrt( 1.0 + temp^2 );
      c = temp * s;
   else
      temp = b / a;
      c = 1.0 / sqrt( 1.0 + temp^2 );
      s = temp * c;
   end

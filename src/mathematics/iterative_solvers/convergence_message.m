% This file is a modification of the file bicgstab.m, downloaded from
% http://it.mathworks.com/matlabcentral/fileexchange/2158-templates-for-the-solution-of-linear-systems
%
% Amos Egel (KIT) made the following changes to the original file: 
% - convergence monitor output to console
% - input matrices A and M are function handles
% - change order of input arguments to match Matlabs built-in bicgstab
% - change check of half step convergency to norm(s)/bnrm2 instead
%   of only norm(s)
%
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

%======================================================================
%> @brief Display convergence status to console
%>
%> @param msg (string): Previous convergence message
%> @param iternum (int): Iteration number
%> @param resid (float): Residual
%>
%> @retval msg (string): convergence message displayed to screen
%======================================================================
function msg = convergence_message(msg, iternum, resid)

    % fprintf(repmat('\b', [1, length(msg)]));
    msg = [sprintf('iteration number: %i, relative residual: %f', iternum, resid)];
    fprintf('%s\n', msg);

end

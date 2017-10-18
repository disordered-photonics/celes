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

%===============================================================================
%> @brief Customized BiCGStab algorithm to allow for live monitoring of
%> convergence progress.
%>
%> This function is a modified version of the algorithm provided with the
%> package "Templates for the Solution of Linear Systems"  
%> from matlab file exchange
%> 
%> @param A (function handle): The linear operator to be solved
%> @param b (float array): Right hand side of linear system
%> @param tol (float): Convergence tolerance - algorithm stops when 
%>        relative error is smaller than tol
%> @param max_it (int): Maximal number of iterations
%> @param Minv (function handle): Preconditioner operator, should be
%>        similar to the inverse of A
%> @param 0 (-): Empty input parameter to make format consistent with
%>        MATLAB's built-in bicgstab
%> @param x (float array): Initial guess for solution vector
%>
%> @retval x (float array): approximate solution of linear system
%> @retval error (float): residual relative error
%> @retval iter (int): number of iterations
%> @retval flag (int): indicates success of algorithm:
%>         0 = solution found to tolerance
%>         1 = no convergence given max_it
%>        -1 = breakdown: rho = 0
%>        -2 = breakdown: omega = 0
%> @retval convHist (float array): relative error as a function of
%>         iteration half-steps
%===============================================================================
function [x, error, iter, flag, convhist] = bicgstab_custom(A, b, tol, max_it, Minv, ~,x)

% original header:
%  -- Iterative template routine --
%     Univ. of Tennessee and Oak Ridge National Laboratory
%     October 1, 1993
%     Details of this algorithm are described in "Templates for the
%     Solution of Linear Systems: Building Blocks for Iterative
%     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra,
%     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
%     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
%
%  [x, error, iter, flag] = bicgstab(A, x, b, M, max_it, tol)
%
% bicgstab.m solves the linear system Ax=b using the
% BiConjugate Gradient Stabilized Method with preconditioning.
%
% input   A        REAL matrix
%         x        REAL initial guess vector
%         b        REAL right hand side vector
%         M        REAL preconditioner matrix
%         max_it   INTEGER maximum number of iterations
%         tol      REAL error tolerance
%
% output  x        REAL solution vector
%         error    REAL error norm
%         iter     INTEGER number of iterations performed
%         flag     INTEGER: 0 = solution found to tolerance
%                           1 = no convergence given max_it
%                          -1 = breakdown: rho = 0
%                          -2 = breakdown: omega = 0
%
%
%

iter = 0;                   % initialization
flag = 0;
convhist = [];
msg='';
bnrm2 = norm( b );
if  ( bnrm2 == 0.0 ), bnrm2 = 1.0; end
r = b - A(x);
error = norm( r ) / bnrm2;
if ( error < tol ), return, end
omega  = 1.0;
r_tld = r;
for iter = 1:max_it
    rho   = r_tld'*r;
    if ( rho == 0.0 ), break, end
    if ( iter > 1 )
        beta  = ( rho/rho_1 )*( alpha1/omega );
        p = r + beta*( p - omega*v );
    else
        p = r;
    end
    p_hat = Minv(p);
    v = A(p_hat);
    alpha1 = rho / ( r_tld'*v );
    s = r - alpha1*v;
    errhalf = norm(s)/bnrm2;

    % status output
    msg=convergence_message(msg,iter-0.5, errhalf);
    convhist(end+1)=gather(errhalf);

    if ( errhalf < tol )                        % early convergence check
        x = x + alpha1*p_hat;
        break;
    end
    s_hat = Minv(s);
    t = A(s_hat);
    omega = ( t'*s) / ( t'*t );

    x = x + alpha1*p_hat + omega*s_hat;         % update approximation

    r = s - omega*t;
    error = norm( r ) / bnrm2;                  % check convergence

    % status output
    msg=convergence_message(msg,iter, errhalf);
    convhist(end+1)=gather(error);

    if ( error <= tol ), break, end
    if ( omega == 0.0 ), break, end
    rho_1 = rho;
end

if ( error <= tol || norm(s)/bnrm2 <= tol )     % converged
    if ( norm(s)/bnrm2 <= tol )
        error = norm(s) / bnrm2;
    end
    flag =  0;
elseif ( omega == 0.0 )                         % breakdown
    flag = -2;
elseif ( rho == 0.0 )
    flag = -1;
else                                            % no convergence
    flag = 1;
end
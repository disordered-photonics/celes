% This file is a modification of the file gmres.m, downloaded from
% http://it.mathworks.com/matlabcentral/fileexchange/2158-templates-for-the-solution-of-linear-systems
%
% Amos Egel (KIT) made the following changes to the original file: 
% - convergence monitor output to console
% - input matrices A and M are function handles
% - change order of input arguments to match Matlabs built-in gmres
% - fix apparent bug in update of x, where x = x + V*y throws an error
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
%> @brief Customized GMRES algorithm to allow for live monitoring of
%> convergence progress.
%>
%> This function is a modified version of the algorithm provided with the
%> package "Templates for the Solution of Linear Systems"  
%> from matlab file exchange
%> 
%> @param A (function handle): The linear operator to be solved
%> @param b (float array): Right hand side of linear system
%> @param restrt (int): Restart GMRES after how many iterations
%> @param tol (float): Convergence tolerance - algorithm stops when 
%> relative error is smaller than tol
%> @param max_it (int): Maximal number of iterations
%> @param Minv (function handle): Preconditioner operator, should be
%> similar to the inverse of A
%> @param ~ (-): Empty input parameter to make format consistent with
%> Matlab's built-in bicgstab
%> @param x (float array): Initial guess for solution vector
%>
%> @retval x (float array): approximate solution of linear system
%> @retval error (float): residual relative error
%> @retval iter (int): number of iterations
%> @retval flag (int): indicates success of algorithm:
%> 0 = solution found to tolerance
%> 1 = no convergence given max_it
%> -1 = breakdown: rho = 0
%> -2 = breakdown: omega = 0
%> @retval convHist (float array): relative error as a function of
%> iteration number
%======================================================================
function [x, error, iter, flag,convHist] = gmres_custom( A, b, restrt, tol, max_it, Minv , ~, x)

% customized GMRES, made by modifying the file "gmres.m" from the package
% Templates for the Solution of Linear Systems
% http://it.mathworks.com/matlabcentral/fileexchange/2158-templates-for-the-solution-of-linear-systems
%
% changes to the original file:
% - convergence monitor output to console
% - input matrices A and M are function handles
% - change order of input arguments to match Matlabs built-in gmres
% - fix apparent bug in update of x, where x = x + V*y throws an error

% original header:
% " -- Iterative template routine --
%     Univ. of Tennessee and Oak Ridge National Laboratory
%     October 1, 1993
%     Details of this algorithm are described in "Templates for the
%     Solution of Linear Systems: Building Blocks for Iterative
%     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra,
%     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
%     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
%
% [x, error, iter, flag] = gmres( A, x, b, M, restrt, max_it, tol )
%
% gmres.m solves the linear system Ax=b
% using the Generalized Minimal residual ( GMRESm ) method with restarts .
%
% input   A        REAL nonsymmetric positive definite matrix
%         x        REAL initial guess vector
%         b        REAL right hand side vector
%         M        REAL preconditioner matrix
%         restrt   INTEGER number of iterations between restarts
%         max_it   INTEGER maximum number of iterations
%         tol      REAL error tolerance
%
% output  x        REAL solution vector
%         error    REAL error norm
%         iter     INTEGER number of iterations performed
%         flag     INTEGER: 0 = solution found to tolerance
%                           1 = no convergence given max_it  "
%
%
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

iter = 0;                                         % initialization
flag = 0;
fprintf(1,'\n');
msg='';
convHist=[];

bnrm2 = norm( b );
if  ( bnrm2 == 0.0 ), bnrm2 = 1.0; end

r = Minv( b-A(x) );
error = norm( r ) / bnrm2;
if ( error < tol ) return, end

[n,~] = size(b);                                  % initialize workspace
m = restrt;
if isa(b,'gpuArray')
    V(1:n,1:m+1) = gpuArray.zeros(n,m+1,'single');
else
    V(1:n,1:m+1) = zeros(n,m+1,'single');
end
H(1:m+1,1:m) = zeros(m+1,m);
cs(1:m) = zeros(m,1);
sn(1:m) = zeros(m,1);
e1    = zeros(n,1);
e1(1) = 1.0;

for iter = 1:max_it                             % begin iteration
    
    r = Minv( b-A(x) );
    V(:,1) = r / norm( r );
    s = gather(norm( r )*e1);
    for i = 1:m                                   % construct orthonormal
        w = Minv(A(V(:,i)));                         % basis using Gram-Schmidt
        for k = 1:i
            H(k,i)= gather(w'*V(:,k));
            w = w - H(k,i)*V(:,k);
        end
        H(i+1,i) = gather(norm( w ));
        V(:,i+1) = w / H(i+1,i);
        for k = 1:i-1                              % apply Givens rotation
            temp     =  cs(k)*H(k,i) + sn(k)*H(k+1,i);
            H(k+1,i) = -sn(k)*H(k,i) + cs(k)*H(k+1,i);
            H(k,i)   = temp;
        end
        [cs(i),sn(i)] = rotmat( H(i,i), H(i+1,i) ); % form i-th rotation matrix
        temp   = cs(i)*s(i);                        % approximate residual norm
        s(i+1) = -sn(i)*s(i);
        s(i)   = temp;
        H(i,i) = cs(i)*H(i,i) + sn(i)*H(i+1,i);
        H(i+1,i) = 0.0;
        error  = abs(s(i+1)) / bnrm2;
        
        fprintf(1,repmat('\b',[1,length(msg)]));
        msg1 = sprintf('iteration %0.1f',(iter-1)*m+i);
        gap1 = repmat(' ',[1,15-length(msg1)]);
        msg2 = sprintf('relative residual %0.2e',error);
        gap2 = repmat(' ',[1,25-length(msg2)]);
        msg = [msg1,gap1,'--',msg2,gap2];
        fprintf(1,msg);
        convHist(end+1)=gather(error);

        if ( error <= tol )                        % update approximation
            y = H(1:i,1:i) \ s(1:i);                 % and exit
            x = x + V(:,1:i)*y;
            break;
        end
    end
    
    if ( error <= tol ), break, end
    y = H(1:m,1:m) \ s(1:m);
    %x = x + V*y;                                   % update approximation
    x = x + V(:,1:m)*y;                                   % update approximation
    r = Minv( b-A(x) );                              % compute residual
    s(i+1) = norm(r);
    error = s(i+1) / bnrm2;                        % check convergence
    if ( error <= tol ), break, end;
end

if ( error > tol ) flag = 1; end;                 % converged


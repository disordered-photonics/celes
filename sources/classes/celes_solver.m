%> @file celes_solver.m
% ======================================================================
%> @brief Solver settings for the master equation 
%> \f$\sum_{i',n'}M^{ii'}_{nn'} b_{n'}^{i'} = \sum_{n'} T^i_{nn'} a_{\mathrm{in},n'}^i\f$
% ======================================================================


classdef celes_solver

    properties
        %> solver type, at the moment 'BiCGStab' and 'GMRES' implemented
        type = 'BiCGStab';
        
        %> solution tolerance. solver stops when relative error smaller 
        tolerance=1e-4;
        
        %> maximal number of iterations. solver stops if exceeded
        maxIter=1000;
        
        %> for type='GMRES': restart parameter
        restart=10;
        
        %> for monitor=true, the progress of the solver can be seen in the
        %> terminal. For that case, a modified version of the solvers
        %> downlowaded from Matlab Stack Exchange will be employed instead
        %> of the matlab built-in algorithms
        monitor=true;
        
        %> celes_preconditioner object for faster convergency of solver
        preconditioner=celes_preconditioner;
    end
    
    properties (Dependent)
    end
    
    methods
        
        % ======================================================================
        %> @brief Run the solver
        %>
        %> @param mmm (function handle): master matrix multiply operator M*
        %> @param rhs (float array): right hand side of linear system
        %> @param Optional: y0 (float array): initial guess for solution
        %> vector
        %> @return y = M\rhs, 
        %> an approximation to the solution of the linear system
        % ======================================================================
        function [value,convergenceHistory] = run(obj,mmm,rhs,varargin)
            prh = @(x) obj.preconditioner.run(x);
            if isempty(varargin)
                initial_guess=rhs;
            else
                initial_guess=varargin{1};
            end
            switch obj.type
                case 'BiCGStab'
                    if obj.monitor
                        [value,~,~,~,convergenceHistory] = bicgstab_custom(mmm,rhs,obj.tolerance,obj.maxIter,prh,[],initial_guess);
                    else
                        [value,~,~,~,convergenceHistory] = bicgstab(mmm,rhs,obj.tolerance,obj.maxIter,prh,[],initial_guess);
                    end
                case 'GMRES'
                    if obj.monitor
                        [value,~,~,~,convergenceHistory] = gmres_custom(mmm,rhs,obj.restart,obj.tolerance,obj.maxIter,prh,[],initial_guess);
                    else
                        [value,~,~,~,convergenceHistory] = gmres(mmm,rhs,obj.restart,obj.tolerance,obj.maxIter,prh,[],initial_guess);
                    end
            end
        end
    end
end


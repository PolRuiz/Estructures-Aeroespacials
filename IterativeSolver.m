classdef IterativeSolver < Solver
    %Iteratively solve equation system using pcg function    
    properties
        %
    end
    
    methods (Static)
        function result = systSolve(LHS,RHS)
            result = pcg(LHS,RHS,[],100);
            disp('Iterative');
        end
    end
end


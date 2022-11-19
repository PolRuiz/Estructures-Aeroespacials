classdef IterativeSolver < Solver
    %Iteratively solve equation system using pcg function
    
    methods (Access = public)
        function obj = IterativeSolver(cParams)
            obj.init(cParams);
        end
        
        function result = solveSystem(obj)
            result = pcg(obj.LHS,obj.RHS,[],100);
            %disp('Iterative');
        end
    end
end


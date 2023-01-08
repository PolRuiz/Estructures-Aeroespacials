classdef DirectSolver < Solver
    methods (Access = public)
        function obj = DirectSolver(cParams)
            obj.init(cParams);
        end        
        function result = solveSystem(obj)
            result = obj.LHS\obj.RHS;
        end
    end
end


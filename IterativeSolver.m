classdef IterativeSolver < Solver    
    methods (Access = public)
        function obj = IterativeSolver(cParams)
            obj.init(cParams);
        end
        
        function result = solveSystem(obj)
            result = pcg(obj.LHS,obj.RHS,[1e-13],100);
        end
    end
end


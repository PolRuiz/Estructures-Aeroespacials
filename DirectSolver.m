classdef DirectSolver < Solver
    %Directly solve equation system
    
    methods (Access = public)
        function obj = DirectSolver(cParams)
            obj.init(cParams);
        end
        
        function result = solveSystem(obj)
            result = obj.LHS\obj.RHS;
            %disp('Direct');
        end
    end
end


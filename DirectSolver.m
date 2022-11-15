classdef DirectSolver < Solver
    %Directly solve equation system
    properties
        %
    end
    
    methods (Static)
        %no object constructor since we are using static methods
        function result = systSolve(LHS,RHS)
            result = LHS\RHS;
            disp('Direct');
        end
    end
end


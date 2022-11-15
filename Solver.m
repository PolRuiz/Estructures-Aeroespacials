classdef Solver < handle
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods (Static)
        function type = selectSolver(solver)
            switch solver
                case 'Direct'
                    type = DirectSolver;
                case 'Iterative'
                    type = IterativeSolver;
            end
        end
    end
end


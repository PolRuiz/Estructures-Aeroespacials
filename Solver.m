classdef Solver < handle
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = protected)
        RHS
        LHS
    end
    
    methods (Access = public)
        function uL = solve(obj)
            uL = obj.solveSystem();
        end
    end
    
    methods (Static)
        function solver = create(cParams)
            switch cParams.type
                case 'Direct'
                    solver = DirectSolver(cParams);
                case 'Iterative'
                    solver = IterativeSolver(cParams);
            end
        end
    end
    
    methods (Access = protected)
        function init(obj,cParams)
            obj.RHS = cParams.RHS;
            obj.LHS = cParams.LHS;
        end
    end
    
    
    methods (Abstract)
        solveSystem();
    end
    
    
end


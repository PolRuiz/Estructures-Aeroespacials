classdef Solver < handle
    
    properties (Access = protected)
        RHS
        LHS
    end
    
    methods (Access = public)
        function uL = solve(obj)
            uL = obj.solveSystem();            
        end
    end
    
    methods (Access = public, Static)
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
            obj.RHS = cParams.system.RHS;
            obj.LHS = cParams.system.LHS;
        end
    end
    
    
    methods (Access = public, Abstract)
        solveSystem();
    end
    
    
end


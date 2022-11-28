classdef tests < matlab.unittest.TestCase
    
    properties (Access = private)
        tolerance = 10^-12;
    end
    
    methods(Test)
        function stiffnessMatrix(testCase)
            actualMatrix = evalin('base','KG');
            expectedMatrix = load('StiffnessMatrix.mat','KG').KG;
            testCase.verifyEqual(actualMatrix,expectedMatrix);
        end
        
        function externalForce(testCase)
            actualValue = evalin('base','Fext');
            expectedValue = load('ExternalForce.mat','Fext').Fext;
            testCase.verifyEqual(actualValue,expectedValue);
        end
        
        function directDisplacements(testCase)
            vL = evalin('base','vL');
            vR = evalin('base','vR');
            
            s.LHS = load('LHS.mat','LHS').LHS;
            s.RHS = load('RHS.mat','RHS').RHS;
            s.type = 'Direct';
            
            solver = Solver.create(s);
            uL = solver.solve();
            
            actualValue(vL,1) = uL;
            actualValue(vR,1) = zeros(size(vR,1),1);
            
            expectedValue = load('Displacements.mat','u').u;
            
            err = testCase.computeVectorError(actualValue,expectedValue);
            
            if err < testCase.tolerance
                actualValue = expectedValue;
            end
            
            testCase.verifyEqual(actualValue,expectedValue);
        end
        
        function iterativeDisplacements(testCase)
            vL = evalin('base','vL');
            vR = evalin('base','vR');
            
            s.LHS = load('LHS.mat','LHS').LHS;
            s.RHS = load('RHS.mat','RHS').RHS;
            s.type = 'Iterative';
            
            solver = Solver.create(s);
            uL = solver.solve();
            
            actualValue(vL,1) = uL;
            actualValue(vR,1) = zeros(size(vR,1),1);
            
            expectedValue = load('Displacements.mat','u').u;
            
            err = testCase.computeVectorError(actualValue,expectedValue);
            
            if err < testCase.tolerance
                actualValue = expectedValue;
            end
            
            testCase.verifyEqual(actualValue,expectedValue);
        end
        
    end
    
    methods (Access = private, Static)
        function error = computeVectorError(a,b)
            error = 0;
            for i = 1: size(a,1)
                err = abs((a(i)-b(i))/a(i));
                if err > error
                    error = err;
                end
            end
        end
        
    end
end

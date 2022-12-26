classdef tests < matlab.unittest.TestCase
    
    properties (Access = private)
        tolerance = 10^-12;
        
    end
    
    methods(Test)
        function stiffnessMatrix(testCase)
            dimensionalData = load('dimensionalData.mat').dimensionalData;
            s.n_el = dimensionalData.n_el;
            s.n_dof = dimensionalData.n_dof;
            s.n_el_dof = dimensionalData.n_el_dof;
            
            s.Tn = load('mesh.mat').mesh.nodalConnec;            
            s.Td = load('dofConnec.mat').DOFConnec;
            s.Kel = load('elementStiffness.mat').elementStiffness;
            
            c = GlobalStiffnessMatrixComputer(s);
            
            actualMatrix = c.compute();
            expectedMatrix = load('stiffnessMatrix.mat').KG;
            testCase.verifyEqual(actualMatrix,expectedMatrix);
        end
        
        function externalForce(testCase)
            dimensionalData = load('dimensionalData.mat').dimensionalData;
            s.n = dimensionalData.n;
            s.n_d = dimensionalData.n_d;
            s.n_el = dimensionalData.n_el;
            s.n_dof = dimensionalData.n_dof;
            
            s.mesh = load('mesh.mat').mesh;
            s.materialData = load('materialData.mat').materialData;
            
            parameters = load('parameters.mat').parameters;
            s.Wm = parameters.hangingMass;
            s.AeroM = parameters.aeroMultiplier;
            
            geometricalData = load('geometricalData.mat').geometricalData;
            
            s.W = geometricalData.W;
            s.H = geometricalData.H;
            s.D1 = geometricalData.D1;
            s.d1 = geometricalData.d1;
            s.D2 = geometricalData.D2;
            
            c = GlobalForceComputer(s);
            actualValue = c.compute();
            expectedValue = load('externalForce.mat').Fext;
            testCase.verifyEqual(actualValue,expectedValue);
        end
        
        function directDisplacements(testCase)
            vL = load('freeDOF.mat').vL;
            vR = load('fixedDOF.mat').vR;
            
            sys.LHS = load('LHS.mat','LHS').LHS;
            sys.RHS = load('RHS.mat','RHS').RHS;
            
            s.type = 'Direct';
            s.system = sys;
            
            solver = Solver.create(s);
            uL = solver.solve();
            
            actualValue(vL,1) = uL;
            actualValue(vR,1) = zeros(size(vR,1),1);
            
            expectedValue = load('displacements.mat','u').u;
            
            err = testCase.computeVectorError(actualValue,expectedValue);
            
            if err < testCase.tolerance
                actualValue = expectedValue;
            end
            
            testCase.verifyEqual(actualValue,expectedValue);
        end
        
        function iterativeDisplacements(testCase)
            vL = load('freeDOF.mat').vL;
            vR = load('fixedDOF.mat').vR;
            
            sys.LHS = load('LHS.mat').LHS;
            sys.RHS = load('RHS.mat').RHS;
            
            s.type = 'Iterative';
            s.system = sys;
            
            solver = Solver.create(s);
            uL = solver.solve();
            
            actualValue(vL,1) = uL;
            actualValue(vR,1) = zeros(size(vR,1),1);
            
            expectedValue = load('displacements.mat').u;
            
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

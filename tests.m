classdef tests < matlab.unittest.TestCase
    
    properties (Access = private)
        tolerance = 10^-12;
        
    end
    
    methods(Test)
        function stiffnessMatrix(testCase)
            d = load('dimensionalData.mat').dimensionalData;
            s.n_el     = d.n_el;
            s.n_dof    = d.n_dof;
            s.n_el_dof = d.n_el_dof;
            
            s.Tn  = load('mesh.mat').mesh.nodalConnec;            
            s.Td  = load('dofConnec.mat').DOFConnec;
            s.Kel = load('elementStiffness.mat').elementStiffness;            
            c = GlobalStiffnessMatrixComputer(s);
            
            actualMatrix = c.compute();
            expectedMatrix = load('stiffnessMatrix.mat').KG;
            testCase.verifyEqual(actualMatrix,expectedMatrix);
        end
        
        function externalForce(testCase)
            d = load('dimensionalData.mat').dimensionalData;
            s.n = d.n;
            s.n_d = d.n_d;
            s.n_el = d.n_el;
            s.n_dof = d.n_dof;
            
            s.mesh         = load('mesh.mat').mesh;
            s.materialData = load('materialData.mat').materialData;
            
            parameters = load('parameters.mat').parameters;
            s.hMass = parameters.hangingMass;
            s.AeroM = parameters.aeroMultiplier;
            
            gD = load('geometricalData.mat').geometricalData;
            
            s.W = gD.W;
            s.H = gD.H;
            s.D1 = gD.D1;
            s.d1 = gD.d1;
            s.D2 = gD.D2;
            
            c = GlobalForceComputer(s);
            actualValue = c.compute();
            expectedValue = load('externalForce.mat').externalForce;
            testCase.verifyEqual(actualValue,expectedValue);
        end
        
        function directDisplacements(testCase)

            
            sys.LHS = load('LHS.mat','LHS').LHS;
            sys.RHS = load('RHS.mat','RHS').RHS;
            
            s.type = 'Direct';
            s.system = sys;
            
            solver = Solver.create(s);
            uL = solver.solve();

            vL = load('freeDOF.mat').vL;
            vR = load('fixedDOF.mat').vR;            
            
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

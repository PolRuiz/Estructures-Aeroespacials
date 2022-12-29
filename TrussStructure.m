classdef TrussStructure < handle

    properties (Access = private)
        parameters
        geometricalData
        mesh
        materialData
        fixNod
        dimensionalData
    end

    properties (Access = private)
        DOFConnec
        elementStiffness

        StiffnessMatrix
        ExternalForces
        
        postprocessData
        reactions
        displacements
        problemDOFManager

    end

    methods (Access = public)

        function obj = TrussStructure(cParams)
            obj.init(cParams);
        end

        function computeProblem(obj)
            obj.computeDOFConnectivity();
            obj.computeElementStiffnes();
            obj.computeStiffnessMatrix();
            obj.computeGlobalForceVector();
            obj.createDOFManager();
            obj.solveSystem();
            obj.computeStrainStress();
            obj.plotPostprocess();
        end
    end

    methods (Access = private)


        function computeDOFConnectivity(obj)
            s.n_el = obj.dimensionalData.n_el;
            s.n_nod = obj.dimensionalData.n_nod;
            s.n_i = obj.dimensionalData.n_i;
            s.Tn = obj.mesh.nodalConnec;
            
            c = DOFConnecter(s);
            obj.DOFConnec = c.compute();
        end

        function computeElementStiffnes(obj)
            s.n_d = obj.dimensionalData.n_d;
            s.n_el = obj.dimensionalData.n_el;
            s.mesh = obj.mesh;
            s.materialData = obj.materialData;
            
            c = ElementStiffnessComputer(s);
            obj.elementStiffness = c.compute();
        end

        function computeStiffnessMatrix(obj)

            s.n_dof = obj.dimensionalData.n_dof;
            s.n_el = obj.dimensionalData.n_el;
            s.n_el_dof = obj.dimensionalData.n_el_dof;

            s.Tn = obj.mesh.nodalConnec;            
            s.Kel = obj.elementStiffness;
            s.Td = obj.DOFConnec;

            c = GlobalStiffnessMatrixComputer(s);
            obj.StiffnessMatrix = c.compute();             
        end

        function computeGlobalForceVector(obj)
            s.mesh = obj.mesh;
            s.materialData = obj.materialData;
            
            gD = obj.geometricalData;
            s.W  = gD.W;
            s.H  = gD.H;
            s.D1 = gD.D1;
            s.d1 = gD.d1;
            s.D2 = gD.D2;            
            s.n     = obj.dimensionalData.n;
            s.n_d   = obj.dimensionalData.n_d;
            s.n_el  = obj.dimensionalData.n_el;
            s.n_dof = obj.dimensionalData.n_dof;
            
            s.Wm    = obj.parameters.hangingMass;
            s.AeroM = obj.parameters.aeroMultiplier;
            
            c = GlobalForceComputer(s);
            obj.ExternalForces = c.compute();
        end

        function createDOFManager(obj)
            s.fixNod = obj.fixNod;
            s.n_dof = obj.dimensionalData.n_dof;
            obj.problemDOFManager = DOFManager(s);
        end

        function solveSystem(obj)
            solverType = 'Direct';
            newMatrices = obj.problemDOFManager.splitMatrix(obj.StiffnessMatrix);
            newVectors = obj.problemDOFManager.splitVector(obj.ExternalForces);

            s.type = solverType;
            s.system = obj.problemDOFManager.constructSystem(newMatrices,newVectors);

            solver = Solver.create(s);

            uL = solver.solve();
            uR = obj.fixNod(:,3);
            RR = newMatrices.matRR*uR+newMatrices.matRL*uL-newVectors.vecR;


            sol.displacements = obj.problemDOFManager.joinVector(uL,uR);
            sol.reactions = obj.problemDOFManager.joinVector(0,RR);

            obj.resultData = sol;
        end

        function computeStrainStress(obj)
            s.n_d   = obj.dimensionalData.n_d;
            s.n_i   = obj.dimensionalData.n_i;
            s.n_nod = obj.dimensionalData.n_nod;
            s.n_el  = obj.dimensionalData.n_el;             
            s.mesh = obj.mesh;
            s.materialData = obj.materialData;            
            s.DOFConnec = obj.DOFConnec;
            s.u = obj.resultData.displacements;
            
            c = StrainStressComputer(s);
            obj.postprocessData = c.compute();
        end

        function plotPostprocess(obj)
            % Plot deformed structure with stress of each bar
            scale = 30; % Adjust this parameter for properly visualizing the deformation
            s.mesh = obj.mesh;
            s.u = obj.resultData.displacements;
            s.sig = obj.postprocessData.stress;            
            plotBarStress3D(s,scale);
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.parameters = cParams.parameters ; % mass hanging from bar
            obj.geometricalData = cParams.geometricalData;
            obj.mesh = cParams.mesh;
            obj.materialData = cParams.materialData;
            obj.fixNod = cParams.fixNod;
            obj.dimensionalData = cParams.dimensionalData;
        end

    end

end
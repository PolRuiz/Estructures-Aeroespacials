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
        stiffnessMatrix
        externalForces
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
            obj.computeDisplacements();
            obj.computeStrainStress();
            obj.plotPostprocess();
        end
    end

    methods (Access = private)
        function init(obj, cParams)
            obj.parameters      = cParams.parameters;
            obj.geometricalData = cParams.geometricalData;
            obj.mesh            = cParams.mesh;
            obj.materialData    = cParams.materialData;
            obj.fixNod          = cParams.fixNod;
            obj.dimensionalData = cParams.dimensionalData;
        end

        function computeDOFConnectivity(obj)
            dD = obj.dimensionalData;
            s.n_el  = dD.n_el;
            s.n_nod = dD.n_nod;
            s.n_i   = dD.n_i;
            s.Tn    = obj.mesh.nodalConnec;
            c = DOFConnecter(s);
            obj.DOFConnec = c.compute();
        end

        function computeElementStiffnes(obj)
            dD  = obj.dimensionalData;
            s.n_d           = dD.n_d;
            s.n_el          = dD.n_el;
            s.n_el_dof      = dD.n_el_dof;
            s.mesh          = obj.mesh;
            s.materialData  = obj.materialData;
            c = ElementStiffnessComputer(s);
            obj.elementStiffness = c.compute();
        end

        function computeStiffnessMatrix(obj)
            dD = obj.dimensionalData;
            s.n_dof     = dD.n_dof;
            s.n_el      = dD.n_el;
            s.n_el_dof  = dD.n_el_dof;
            s.Tn        = obj.mesh.nodalConnec;
            s.Kel       = obj.elementStiffness;
            s.Td        = obj.DOFConnec;
            c = GlobalStiffnessMatrixComputer(s);
            obj.stiffnessMatrix = c.compute();
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
            s.hMass = obj.parameters.hangingMass;
            s.AeroM = obj.parameters.aeroMultiplier;
            c = GlobalForceComputer(s);
            obj.externalForces = c.compute();
        end

        function createDOFManager(obj)
            s.fixNod = obj.fixNod;
            s.n_dof = obj.dimensionalData.n_dof;
            obj.problemDOFManager = DOFManager(s);
        end

        function computeDisplacements(obj)
            s.KG        = obj.stiffnessMatrix;
            s.extF      = obj.externalForces;
            s.manager   = obj.problemDOFManager;
            s.uR        = obj.fixNod(:,3);
            c = DisplacementsComputer(s);
            obj.displacements = c.compute();
        end

        function computeStrainStress(obj)
            dD = obj.dimensionalData;
            s.n_d       = dD.n_d;
            s.n_i       = dD.n_i;
            s.n_nod     = dD.n_nod;
            s.n_el      = dD.n_el;
            s.n_el_dof  = dD.n_el_dof;
            s.mesh      = obj.mesh;
            s.Td        = obj.DOFConnec;
            s.mData     = obj.materialData;
            s.u         = obj.displacements;
            c = StrainStressComputer(s);
            obj.postprocessData = c.compute();
        end

        function plotPostprocess(obj)
            scale   = 30;
            s.mesh  = obj.mesh;
            s.u     = obj.displacements;
            s.sig   = obj.postprocessData.stress;
            plotBarStress3D(s,scale);
        end

    end
end
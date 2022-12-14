classdef TrussStructure < handle

    properties (Access = public)
        KG
        Fext
        vectorialData
        preprocessData
        resultData
    end

    properties (Access = private)
        hangingMass
        aeroMultiplier
        geometricalData
        dimensionalData
        DOFConnec
        elementStiffness 
    end


    methods (Access = public)

        function obj = TrussStructure(cParams)
            obj.init(cParams);
        end

        function createData(obj)
            obj.createGeometry();
            obj.computePreprocess();
            obj.createDimensions();
        end

        function computeDOFConnectivity(obj)
            s.preprocessData = obj.preprocessData; %no utilitzo tot
            s.dimensionalData = obj.dimensionalData; %no utilitzo tot
            c = DOFConnecter(s);
            obj.DOFConnec = c.compute();
        end

        function computeElementStiffnes(obj)
            s.preprocessData = obj.preprocessData; %no utilitzo tot
            s.dimensionalData = obj.dimensionalData; %no utilitzo tot
            c = ElementStiffnessComputer(s);
            obj.elementStiffness = c.compute();
        end

        function computeStiffnessMatrix(obj)
            s.preprocessData = obj.preprocessData; %no utilitzo tot
            s.dimensionalData = obj.dimensionalData; %no utilitzo tot
            s.Kel = obj.elementStiffness;
            s.Td = obj.DOFConnec;
            c = GlobalStiffnessMatrixComputer(s);
            obj.KG = c.compute();
        end

        function computeGlobalForceVector(obj)
            s.preprocessData = obj.preprocessData;
            s.geometricalData = obj.geometricalData;
            s.dimensionalData = obj.dimensionalData;
            s.hangingMass = obj.hangingMass;
            s.aeroMultiplier = obj.aeroMultiplier;
            c = GlobalForceComputer(s);
            obj.Fext = c.compute();
        end

        function applyBoundaryCond(obj)
            fixNod = obj.preprocessData.fixNode;
            n_dof = obj.dimensionalData.n_dof;

            c=size(fixNod,1);
            s.uR=fixNod(:,3);
            s.vR=fixNod(:,2);
            DOF=transpose(1:1:n_dof);
            for k=1:n_dof
                for j=1:c
                    if DOF(k)==fixNod(j,2)
                        DOF(k)=0;
                    end
                end
            end
            s.vL=nonzeros(DOF);
            obj.vectorialData = s;
        end

        function computeStrainStress(obj, displacements, freeDisplacements)
            obj.vectorialData.u = displacements;
            obj.vectorialData.uL = freeDisplacements;

            s.dimensionalData = obj.dimensionalData;
            s.vectorialData = obj.vectorialData;
            s.preprocessData = obj.preprocessData;
            s.DOFConnec = obj.DOFConnec;

            c = StrainStressComputer(s);
            obj.resultData = c.compute();
        end
    end

    methods (Access = private)

        function init(obj, cParams)
            obj.hangingMass = cParams.Mass ; % mass hanging from bar
            obj.aeroMultiplier = cParams.AeroM;
        end

        function createGeometry(obj)
            geo.H = 0.9;
            geo.W = 0.85;
            geo.B = 3.2;
            geo.D1 = 18/1000;
            geo.d1 = 7.5/1000;
            geo.D2 = 3/1000;
            obj.geometricalData = geo;
        end

        function computePreprocess(obj)
            s.geometricalData = obj.geometricalData;
            c = PreprocessComputer(s);
            obj.preprocessData = c.compute();
        end

        function createDimensions(obj)
            dim.n_d = size(obj.preprocessData.coor,2);              % Number of dimensions
            dim.n_i = dim.n_d;                                          % Number of DOFs for each node
            dim.n = size(obj.preprocessData.coor,1);                % Total number of nodes
            dim.n_dof = dim.n_i*dim.n;                                      % Total number of degrees of freedom
            dim.n_el = size(obj.preprocessData.nodalConnec,1);      % Total number of elements
            dim.n_nod = size(obj.preprocessData.nodalConnec,2);     % Number of nodes for each element
            dim.n_el_dof = dim.n_i*dim.n_nod;                               % Number of DOFs for each element
            obj.dimensionalData = dim;
        end

%         function solveSystem(obj)
%             %DOFManager
%         end




    end





end
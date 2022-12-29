classdef DataCreator < handle
    
    properties (Access = private)
        geometricalData
        mesh
        materialData
        fixNod
        dimensionalData
        parameters
    end
    
    methods (Access = public)
        
        function problemData = createData(obj)
            obj.inputProblemParameters();
            obj.createGeometry();
            obj.computeMesh();
            obj.computeMaterial();
            obj.createFixNodes();
            obj.createDimensions();
            
            s.parameters = obj.parameters;
            s.geometricalData = obj.geometricalData;
            s.mesh = obj.mesh;
            s.materialData = obj.materialData;
            s.fixNod = obj.fixNod;
            s.dimensionalData = obj.dimensionalData;
            problemData = s;
        end
    end
    
    methods (Access = private)
        function inputProblemParameters(obj)
            s.hangingMass = 150; %kg
            s.aeroMultiplier = 1;
            obj.parameters = s;
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

        function computeMesh(obj)
            s.geometricalData = obj.geometricalData;
            c = MeshComputer(s);
            obj.mesh = c.compute();            
        end
        
        function computeMaterial(obj)
            s.geometricalData = obj.geometricalData;
            c = MaterialComputer(s);
            obj.materialData = c.compute();
        end
        
        function createFixNodes(obj)
            obj.fixNod = [4 10 0; 4 11 0; 4 12 0; 5 13 0; 5 15 0; 7 21 0];
        end
        
        function createDimensions(obj)
            dim.n_d = size(obj.mesh.coor,2);              % Number of dimensions
            dim.n_i = dim.n_d;                                          % Number of DOFs for each node
            dim.n = size(obj.mesh.coor,1);                % Total number of nodes
            dim.n_dof = dim.n_i*dim.n;                                      % Total number of degrees of freedom
            dim.n_el = size(obj.mesh.nodalConnec,1);      % Total number of elements
            dim.n_nod = size(obj.mesh.nodalConnec,2);     % Number of nodes for each element
            dim.n_el_dof = dim.n_i*dim.n_nod;                               % Number of DOFs for each element
            obj.dimensionalData = dim;
        end
        
        
    end
end



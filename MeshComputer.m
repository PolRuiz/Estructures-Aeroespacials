classdef MeshComputer < handle
    
    
    properties (Access = private)
        geometricalData
        mesh
    end
    
    methods (Access = public)
        function obj = MeshComputer(cParams)
            obj.init(cParams);
        end
        
        function mesh = compute(obj)
            obj.computeNodalCoordinates();
            obj.createNodalConnectivities();
            mesh = obj.mesh;
        end
    end
    
    
    methods (Access = private)
        function init(obj, cParams)
            obj.geometricalData = cParams.geometricalData;
        end
        
        function computeNodalCoordinates(obj)
            H = obj.geometricalData.H;
            W = obj.geometricalData.W;
            B = obj.geometricalData.B;
            
            x = [
                2*W,  -W/2,     0;
                2*W,   W/2,     0;
                2*W,     0,     H;
                0,     0,     H;
                0,    -B,     H;
                0,     B,     H;
                W,     0,     H;
                ];
            obj.mesh.coor = x;
        end
        
        function createNodalConnectivities(obj)
            Tn = [1 2; 2 3; 1 3; 3 5; 3 6; 4 6; 6 7; 4 5;
                5 7; 3 7; 4 7; 1 5; 2 6; 1 7; 2 7; 1 4; 2 4];
            obj.mesh.nodalConnec = Tn;
        end
        
    end
end


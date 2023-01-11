classdef MaterialComputer < handle

   
    properties (Access = private)
        materialData
        geometricalData
    end
    
    methods (Access = public)
        function obj = MaterialComputer(cParams)
            obj.init(cParams);            
        end
        
        function materialData = compute(obj)
            obj.computeMaterialProperties();
            obj.createMaterialConnectivities();
            materialData = obj.materialData;
        end
    end
    
    methods (Access = private)
        function init(obj, cParams)
            obj.geometricalData = cParams.geometricalData;
        end
        
        function computeMaterialProperties(obj)
            d1 = obj.geometricalData.d1;
            D1 = obj.geometricalData.D1;
            D2 = obj.geometricalData.D2;
            
            mat = [
                75000000000,    pi*(D1/2)^2-pi*(d1/2)^2,    3350,      pi/32*(D1^4-d1^4);%mat(1)
                147000000000,                pi*(D2/2)^2,     950,      pi/32*D2^4       %mat(2)
                ];
            obj.materialData.matProp = mat;
        end
        
        function createMaterialConnectivities(obj)
            Tmat = [1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2];
            obj.materialData.matConnec = Tmat;
        end
    end
end


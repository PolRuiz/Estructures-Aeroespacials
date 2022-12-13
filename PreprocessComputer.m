classdef PreprocessComputer < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
   
    properties (Access = public)
        preprocessData
    end
   
    properties (Access = private)
        geometricalData
    end
   
   
   
    methods (Access = public)
        function obj = PreprocessComputer(cParams)
            obj.init(cParams);
        end
       
        function preprocessData = compute(obj)
            obj.nodalCoordinatesCreator();
            obj.nodalConnectivitiesCreator();
            obj.matPropComputer();
            obj.matConnectivitiesCreator();
            obj.fixNodeCreator();
            preprocessData = obj.preprocessData;
        end
    end
   
    methods (Access = private)
        function init(obj, cParams)
            obj.geometricalData = cParams.geometricalData;
        end
       
        function nodalCoordinatesCreator(obj)
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
            obj.preprocessData.coor = x;
        end
       
        function nodalConnectivitiesCreator(obj)
            Tn = [1 2; 2 3; 1 3; 3 5; 3 6; 4 6; 6 7; 4 5;
                5 7; 3 7; 4 7; 1 5; 2 6; 1 7; 2 7; 1 4; 2 4];
            obj.preprocessData.nodalConnec = Tn;
        end
       
        function matPropComputer(obj)
            d1 = obj.geometricalData.d1;
            D1 = obj.geometricalData.D1;
            D2 = obj.geometricalData.D2;
           
            mat = [
                %Young M.        Section A.                  Density    Inertia
                75000000000,    pi*(D1/2)^2-pi*(d1/2)^2,    3350,      pi/32*(D1^4-d1^4);%mat(1)
                147000000000,                pi*(D2/2)^2,     950,      pi/32*D2^4       %mat(2)
                ];
            obj.preprocessData.matProp = mat;
        end
       
        function matConnectivitiesCreator(obj)
            Tmat = [1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2];
            obj.preprocessData.matConnec = Tmat;
        end
       
        function fixNodeCreator(obj)
            %  fixNod(k,1) = node at which some DOF is prescribed
            %  fixNod(k,2) = DOF prescribed
            %  fixNod(k,3) = prescribed displacement in the corresponding DOF (0 for fixed)
            fixNod = [4 10 0; 4 11 0; 4 12 0; 5 13 0; 5 15 0; 7 21 0];
            obj.preprocessData.fixNode = fixNod;
        end      
       
    end
   
   
end

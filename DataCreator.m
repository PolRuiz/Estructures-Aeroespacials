classdef DataCreator < handle

    properties (Access = public)
        geometricalData
        preprocessData
    end

    methods (Access = public)

        function problemData = createData(obj)
            obj.createGeometry();
            obj.computePreprocess();
            s.geometricalData = obj.geometricalData;
            s.preprocessData = obj.preprocessData;
            problemData = s;
        end
    end

    methods (Access = private)
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


    end
end



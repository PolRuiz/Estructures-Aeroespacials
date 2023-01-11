classdef DisplacementsComputer < handle
    properties(Access = private)
        KG
        extF
        manager
        uR
    end

    methods(Access = public)
        function obj = DisplacementsComputer(cParams)
            obj.init(cParams);
        end

        function displacements = compute(obj)
            solverType  = 'Direct';
            newMatrices = obj.manager.splitMatrix(obj.KG);
            newVectors  = obj.manager.splitVector(obj.extF);
            s.type      = solverType;
            s.system    = obj.manager.constructSystem(newMatrices,newVectors);
            solver      = Solver.create(s);
            uL          = solver.solve();
            displacements   = obj.manager.joinVector(uL,obj.uR);
        end
    end

    methods (Access = private)
        function init(obj, cParams)
            obj.KG      = cParams.KG;
            obj.extF    = cParams.extF;
            obj.manager = cParams.manager;
            obj.uR      = cParams.uR;            
        end

    end
end
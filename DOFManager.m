classdef DOFManager < handle

    properties (Access = private)
        fixNod
        n_dof        
        vR
        vL
        uR
    end

    methods (Access = public)

        function obj = DOFManager(cParams)
            obj.init(cParams);
            obj.createVectors();
        end

        function newMatrices = splitMatrix(obj,matrix)
            s.matLL = matrix(obj.vL,obj.vL);
            s.matLR = matrix(obj.vL,obj.vR);
            s.matRL = matrix(obj.vR,obj.vL);
            s.matRR = matrix(obj.vR,obj.vR);
            newMatrices = s;
        end

        function newVector = splitVector(obj,vector)
            s.vecL = vector(obj.vL,1);
            s.vecR = vector(obj.vR,1);
            newVector = s;
        end

        function newVector = joinVector(obj, freeResult, fixedResult)
            newVector(obj.vL,1) = freeResult;
            newVector(obj.vR,1) = fixedResult;
        end

        function system = constructSystem(obj,stiffnessStruct,forceStruct)
            s.LHS = stiffnessStruct.matLL;
            s.RHS = forceStruct.vecL-stiffnessStruct.matLR*obj.uR;
            system = s;
        end
    end


    methods (Access = private)

        function init(obj, cParams)
            obj.fixNod = cParams.fixNod;
            obj.n_dof = cParams.n_dof;
            
        end

        function createVectors(obj)
            c = size(obj.fixNod,1);
            obj.uR = obj.fixNod(:,3);
            obj.vR = obj.fixNod(:,2);
            DOF = transpose(1:1:obj.n_dof);
            for k = 1:obj.n_dof
                for j = 1:c
                    if DOF(k) == obj.fixNod(j,2)
                        DOF(k) = 0;
                    end
                end
            end
            obj.vL = nonzeros(DOF);

        end

    end
end


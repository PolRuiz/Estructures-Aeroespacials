classdef GlobalStiffnessMatrixComputer < handle

    properties (Access = private)
        KG        
        Kel
        Td
        Tn
        n_dof
        n_el
        n_el_dof
    end

    methods (Access = public)
        function obj = GlobalStiffnessMatrixComputer(cParams)
            obj.init(cParams);
        end

        function kG = compute(obj)
            obj.computeKG();
            kG = obj.KG;
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.Kel         = cParams.Kel;
            obj.Td          = cParams.Td;
            obj.Tn          = cParams.Tn;
            obj.n_dof       = cParams.n_dof;
            obj.n_el        = cParams.n_el;
            obj.n_el_dof    = cParams.n_el_dof;
        end


        function computeKG(obj)
            obj.KG = zeros(obj.n_dof,obj.n_dof);
            for e = 1:obj.n_el
                for i = 1:obj.n_el_dof
                    I = obj.Td(e,i);
                    for j       = 1:obj.n_el_dof
                        J       = obj.Td(e,j);
                        obj.KG(I,J) = obj.KG(I,J) + obj.Kel(i,j,e);
                    end
                end
            end
        end
    end
end


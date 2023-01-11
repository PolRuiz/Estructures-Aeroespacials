classdef DOFConnecter < handle


    properties (Access = private)
        n_el
        n_nod
        n_i
        Tn
    end

    methods (Access = public)

        function obj = DOFConnecter(cParams)
            obj.init(cParams)
        end

        function DOFConnec = compute(obj)
            Td = zeros(obj.n_el, obj.n_nod*obj.n_i);
            for i=1:obj.n_el
                for j = 1:obj.n_nod
                    Td(i,3*j-2) =   3*obj.Tn(i,j)-2;
                    Td(i,3*j-1) =   3*obj.Tn(i,j)-1;
                    Td(i,3*j)   =   3*obj.Tn(i,j);
                end
            end
            DOFConnec = Td;
        end

    end



    methods (Access = private)

        function init(obj, cParams)
            obj.n_el    = cParams.n_el;
            obj.n_nod   = cParams.n_nod;
            obj.n_i     = cParams.n_i;
            obj.Tn      = cParams.Tn;
        end

    end


end

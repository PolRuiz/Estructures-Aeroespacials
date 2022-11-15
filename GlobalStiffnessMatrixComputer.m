classdef GlobalStiffnessMatrixComputer < handle
    %   This class is used to compute the Global Stiffnes Matrix of a system
    %   The arguments are going to be the Kel, Tn, n (number of nodes) and n_i (DOF per node)
    
    
    properties
        KG
        Kel
        Tn
        n
        n_i
        n_el
        n_nod
        n_el_dof
        n_dof
        Td
    end
    
    methods (Access = public)
        function obj = GlobalStiffnessMatrixComputer(Kel,Tn,n,n_i)
            obj.init(Kel,Tn,n,n_i);
            obj.compute_n();
            obj.computeTd();
        end
        
        function kG = obtainStiffnessMatrix(obj)
            obj.computeKG();
            kG = obj.KG;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,Kel,Tn,n,n_i)
            obj.Kel = Kel;
            obj.Tn  = Tn;
            obj.n   = n;
            obj.n_i = n_i;
        end
        
        function compute_n(obj)
            obj.n_el     = size(obj.Tn,1);
            obj.n_nod    = size(obj.Tn,2);
            obj.n_el_dof = obj.n_i*obj.n_nod;
            obj.n_dof    = obj.n_i*obj.n;
        end
        
        function computeTd(obj)
            obj.Td  = zeros(obj.n_el, obj.n_nod*obj.n_i);
            for i   = 1:obj.n_el
                for j = 1:obj.n_nod
                    obj.Td(i,3*j-2) =   3*obj.Tn(i,j)-2;
                    obj.Td(i,3*j-1) =   3*obj.Tn(i,j)-1;
                    obj.Td(i,3*j)   =   3*obj.Tn(i,j);
                end
            end
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


classdef ElementStiffnessComputer < handle
    properties (Access = private)
        n_d
        n_el
        n_el_dof
        x
        Tn
        mat
        Tmat
    end
    
    properties (Access = private)
        elCoord
        l_el
        Kelem
    end
    
    methods (Access = public)
        
        function obj = ElementStiffnessComputer(cParams)
            obj.init(cParams)
        end
        
        function elementStiffness = compute(obj)
            Kel = zeros(obj.n_el_dof, obj.n_el_dof,obj.n_el);
            for iElem = 1:obj.n_el
                obj.computeCoord(iElem);
                obj.computeLength();
                obj.computeElementalStiffnessMatrix(iElem);
                Kel(:,:,iElem) = obj.Kelem;
            end
            elementStiffness = Kel;
        end
    end
    
    methods (Access = private)
        
        function init(obj, cParams)
            obj.n_d         = cParams.n_d;
            obj.n_el        = cParams.n_el;
            obj.n_el_dof    = cParams.n_el_dof;
            obj.x           = cParams.mesh.coor;
            obj.Tn          = cParams.mesh.nodalConnec;
            obj.mat         = cParams.materialData.matProp;
            obj.Tmat        = cParams.materialData.matConnec;
        end
        
        function computeCoord(obj, iElem)
            s.x1 = obj.x(obj.Tn(iElem,1),1);
            s.y1 = obj.x(obj.Tn(iElem,1),2);
            s.z1 = obj.x(obj.Tn(iElem,1),3);
            s.x2 = obj.x(obj.Tn(iElem,2),1);
            s.y2 = obj.x(obj.Tn(iElem,2),2);
            s.z2 = obj.x(obj.Tn(iElem,2),3);
            obj.elCoord = s;
        end
        
        function computeLength(obj)
            x1 = obj.elCoord.x1;
            y1 = obj.elCoord.y1;
            z1 = obj.elCoord.z1;
            x2 = obj.elCoord.x2;
            y2 = obj.elCoord.y2;
            z2 = obj.elCoord.z2;            
            obj.l_el = sqrt((x2-x1)^2 + (y2-y1)^2+(z2-z1)^2);
        end
        
        function computeElementalStiffnessMatrix(obj, iElem)
            E = obj.mat(obj.Tmat(iElem),1);
            A = obj.mat(obj.Tmat(iElem),2);
            Kep = (A*E)/obj.l_el*[1 -1;-1 1];
            Re  = obj.computeRotationMatrix();
            obj.Kelem  = transpose(Re)*Kep*Re;
        end
        
        function Re = computeRotationMatrix(obj)
            x1 = obj.elCoord.x1;
            y1 = obj.elCoord.y1;
            z1 = obj.elCoord.z1;
            x2 = obj.elCoord.x2;
            y2 = obj.elCoord.y2;
            z2 = obj.elCoord.z2;
            Re = 1/obj.l_el*[x2-x1 y2-y1 z2-z1 0 0 0;
                0 0 0 x2-x1 y2-y1 z2-z1];
        end
    end
end


classdef StrainStressComputer < handle


    properties (Access = private)
        n_d
        n_i
        n_nod
        n_el
        n_el_dof
        u
        Td
        x
        Tn
        mat
        Tmat
    end

    properties (Access = private)
        elCoord
        l_el
        epsElem
        sigElem
    end

    methods (Access = public)

        function obj = StrainStressComputer(cParams)
            obj.init(cParams)
        end

        function postprocessData = compute(obj)
            eps = zeros(obj.n_el,1);
            sig = zeros(obj.n_el,1);
            for iElem = 1:obj.n_el              
                obj.computeCoord(iElem)
                obj.computeLength();
                obj.computeElementStainStress(iElem);
                eps(iElem) = obj.epsElem;
                sig(iElem) = obj.sigElem;                
            end            
            postprocessData.strain = eps;
            postprocessData.stress = sig;
        end
    end



    methods (Access = private)

        function init(obj, cParams)
            obj.n_d         = cParams.n_d;
            obj.n_i         = cParams.n_i;
            obj.n_nod       = cParams.n_nod;
            obj.n_el        = cParams.n_el;
            obj.n_el_dof    = cParams.n_el_dof;
            obj.u           = cParams.u;
            obj.Td          = cParams.Td;
            obj.x           = cParams.mesh.coor;
            obj.Tn          = cParams.mesh.nodalConnec;
            obj.mat         = cParams.mData.matProp;
            obj.Tmat        = cParams.mData.matConnec;
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

        function computeElementStainStress(obj, iElem)
            E   = obj.mat(obj.Tmat(iElem),1);
            Re  = obj.computeRotationMatrix;
            u_e = zeros(obj.n_el_dof,1);
            for i = 1:obj.n_el_dof
                I = obj.Td(iElem,i);
                u_e(i) = obj.u(I);
            end
            u_p             =   Re*u_e;
            obj.epsElem     =   (1/obj.l_el)*[-1 1]*u_p;
            obj.sigElem     =   E*obj.epsElem;
        end

    end


end
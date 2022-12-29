classdef ElementStiffnessComputer < handle


    properties (Access = private)
        n_d
        n_el
        x
        Tn
        mat
        Tmat
    end

    methods (Access = public)

        function obj = ElementStiffnessComputer(cParams)
            obj.init(cParams)
        end

        function elementStiffness = compute(obj)
            n_nod = size(obj.Tn,2);
            n_el_dof = obj.n_d*n_nod;

            x1 = zeros(obj.n_el,1);
            y1 = zeros(obj.n_el,2);
            z1 = zeros(obj.n_el,3);

            x2 = zeros(obj.n_el,1);
            y2 = zeros(obj.n_el,2);
            z2 = zeros(obj.n_el,3);


%             Re=zeros(n_nod,n_el_dof);
%             Kep=zeros(n_nod,n_nod);
%             Ke=zeros(n_el_dof, n_el_dof);

            Kel = zeros(n_el_dof, n_el_dof,obj.n_el);
            l = zeros(obj.n_el,1);

            for iElem=1:obj.n_el

                obj.computeCoord(iElem)
                obj.computeLength()                
                obj.computeElementalStiffnessMatrix()
                obj.rotateElementalStiffnessMatrix()
                Ke(:,:,iElem) = obj.Kelem;
            end

                E = obj.mat(obj.Tmat(iElem),1);
                A = obj.mat(obj.Tmat(iElem),2);
                x1(iElem) = obj.x(obj.Tn(iElem,1),1);
                y1(iElem) = obj.x(obj.Tn(iElem,1),2);
                z1(iElem) = obj.x(obj.Tn(iElem,1),3);
                x2(iElem) = obj.x(obj.Tn(iElem,2),1);
                y2(iElem) = obj.x(obj.Tn(iElem,2),2);
                z2(iElem) = obj.x(obj.Tn(iElem,2),3);


                l(iElem) = sqrt((x2(iElem)-x1(iElem))^2 + (y2(iElem)-y1(iElem))^2+(z2(iElem)-z1(iElem))^2);

                Re = 1/l(iElem)*[x2(iElem)-x1(iElem) y2(iElem)-y1(iElem) z2(iElem)-z1(iElem) 0 0 0;
                    0 0 0 x2(iElem)-x1(iElem) y2(iElem)-y1(iElem) z2(iElem)-z1(iElem)];


                Kep = (A*E)/l(iElem)*[1 -1;-1 1];

                Ke = transpose(Re)*Kep*Re;

                Kel(:,:,iElem) = Ke(:,:);
            end

            elementStiffness = Kel;

        end

    end



    methods (Access = private)

        function init(obj, cParams)
            obj.n_d = cParams.n_d;
            obj.n_el = cParams.n_el;
            obj.x = cParams.mesh.coor;
            obj.Tn = cParams.mesh.nodalConnec;
            obj.mat = cParams.materialData.matProp;
            obj.Tmat = cParams.materialData.matConnec;
        end

    end


end

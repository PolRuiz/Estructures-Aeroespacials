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

            for e=1:obj.n_el
                E = obj.mat(obj.Tmat(e),1);
                A = obj.mat(obj.Tmat(e),2);
                x1(e) = obj.x(obj.Tn(e,1),1);
                y1(e) = obj.x(obj.Tn(e,1),2);
                z1(e) = obj.x(obj.Tn(e,1),3);
                x2(e) = obj.x(obj.Tn(e,2),1);
                y2(e) = obj.x(obj.Tn(e,2),2);
                z2(e) = obj.x(obj.Tn(e,2),3);


                l(e) = sqrt((x2(e)-x1(e))^2 + (y2(e)-y1(e))^2+(z2(e)-z1(e))^2);

                Re = 1/l(e)*[x2(e)-x1(e) y2(e)-y1(e) z2(e)-z1(e) 0 0 0;
                    0 0 0 x2(e)-x1(e) y2(e)-y1(e) z2(e)-z1(e)];


                Kep = (A*E)/l(e)*[1 -1;-1 1];

                Ke = transpose(Re)*Kep*Re;

                Kel(:,:,e) = Ke(:,:);
            end

            elementStiffness = Kel;

        end

    end



    methods (Access = private)

        function init(obj, cParams)
            obj.n_d = cParams.dimensionalData.n_d;
            obj.n_el = cParams.dimensionalData.n_el;
            obj.x = cParams.preprocessData.coor;
            obj.Tn = cParams.preprocessData.nodalConnec;
            obj.mat = cParams.preprocessData.matProp;
            obj.Tmat = cParams.preprocessData.matConnec;
        end

    end


end

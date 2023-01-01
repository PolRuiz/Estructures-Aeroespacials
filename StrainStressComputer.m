classdef StrainStressComputer < handle


    properties (Access = private)
        n_d
        n_i
        n_nod
        n_el

        u
        Td

        x
        Tn
        mat
        Tmat

    end

    methods (Access = public)

        function obj = StrainStressComputer(cParams)
            obj.init(cParams)
        end

        function postprocessData = compute(obj)
            c=size(obj.Td,2);
            n_el_dof=obj.n_d*obj.n_nod;

            x1=zeros(obj.n_el,1);
            y1=zeros(obj.n_el,2);
            z1=zeros(obj.n_el,3);

            x2=zeros(obj.n_el,1);
            y2=zeros(obj.n_el,2);
            z2=zeros(obj.n_el,3);

            l=zeros(obj.n_el,1);

            Re=zeros(obj.n_nod,n_el_dof,obj.n_el);
            u_e=zeros(c,obj.n_el);
            u_p=zeros(obj.n_nod,obj.n_el);
            eps=zeros(obj.n_el,1);
            sig=zeros(obj.n_el,1);



            for e=1:obj.n_el
                E=obj.mat(obj.Tmat(e),1);
                A=obj.mat(obj.Tmat(e),2);
                x1(e) = obj.x(obj.Tn(e,1),1);
                y1(e) = obj.x(obj.Tn(e,1),2);
                z1(e) = obj.x(obj.Tn(e,1),3);
                x2(e) = obj.x(obj.Tn(e,2),1);
                y2(e) = obj.x(obj.Tn(e,2),2);
                z2(e) = obj.x(obj.Tn(e,2),3);


                l(e)=sqrt((x2(e)-x1(e))^2 + (y2(e)-y1(e))^2+(z2(e)-z1(e))^2);

                Re(:,:,e)=1/l(e)*[x2(e)-x1(e) y2(e)-y1(e) z2(e)-z1(e) 0 0 0;
                    0 0 0 x2(e)-x1(e) y2(e)-y1(e) z2(e)-z1(e)];





                for i=1:obj.n_i*obj.n_nod

                    I=obj.Td(e,i);
                    u_e(i,e)=obj.u(I);

                end

                u_p(:,e)=Re(:,:,e)*u_e(:,e);

                eps(e,1)=(1/l(e))*[-1 1]*u_p(:,e);

                sig(e,1) = E*eps(e,1);
            end

            postprocessData.strain = eps;
            postprocessData.stress = sig;
        end
    end



    methods (Access = private)

        function init(obj, cParams)
            obj.n_d = cParams.n_d;
            obj.n_i = cParams.n_i;
            obj.n_nod = cParams.n_nod;
            obj.n_el = cParams.n_el;

            obj.u = cParams.u;
            obj.Td = cParams.Td;

            obj.x = cParams.mesh.coor;
            obj.Tn = cParams.mesh.nodalConnec;
            obj.mat = cParams.mData.matProp;
            obj.Tmat = cParams.mData.matConnec;

        end

    end


end

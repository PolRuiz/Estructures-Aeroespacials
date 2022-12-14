classdef GlobalForceComputer < handle


    properties (Access = private)
        n
        n_d
        n_el
        n_dof

        x
        Tn
        mat
        Tmat
        Wm
        AeroM

        W
        H
        D1
        d1
        D2
    end

    methods (Access = public)

        function obj = GlobalForceComputer(cParams)
            obj.init(cParams)
        end

        function Fext = compute(obj)

            F_ext=zeros(obj.n_dof,1);


            Pes=zeros(obj.n_dof,1);
            %             massa=zeros(obj.n_dof,1);

            x1=zeros(obj.n_el,1);
            y1=zeros(obj.n_el,2);
            z1=zeros(obj.n_el,3);

            x2=zeros(obj.n_el,1);
            y2=zeros(obj.n_el,2);
            z2=zeros(obj.n_el,3);
            Ptotal=0;

            l = zeros(obj.n_el,1);

            for e=1:obj.n_el
                nodeA = obj.Tn(e,1);
                nodeB = obj.Tn(e,2);
                
                x1(e) = obj.x(nodeA,1);
                y1(e) = obj.x(nodeA,2);
                z1(e) = obj.x(nodeA,3);
                x2(e) = obj.x(nodeB,1);
                y2(e) = obj.x(nodeB,2);
                z2(e) = obj.x(nodeB,3);


                l(e)=sqrt((x2(e)-x1(e))^2 + (y2(e)-y1(e))^2+(z2(e)-z1(e))^2);
                n1=obj.Tn(e,1);
                n2=obj.Tn(e,2);

                if obj.Tmat(e) == 1

                    V=pi*(obj.D1/2)^2*l(e)-pi*(obj.d1/2)^2*l(e);

                else

                    V=(obj.D2/2)^2*pi*l(e);

                end

                Ptotal=Ptotal+V*obj.mat(obj.Tmat(e),3)*9.81;

                Pes(3*n1)=Pes(3*n1)-V*obj.mat(obj.Tmat(e),3)*9.81/2;
                Pes(3*n2)=Pes(3*n2)-V*obj.mat(obj.Tmat(e),3)*9.81/2;

            end
            %
            % massa=-Pes/9.81;
            mTotal=Ptotal/9.81;


            M_PES_7=0;

            for i=1:obj.n
                M_PES_7=M_PES_7+Pes(3*i)*[obj.x(7,1)-obj.x(i,1)];
            end


            L=obj.Wm+Ptotal;
            T=(obj.Wm*obj.W+2/5*L*obj.W+M_PES_7)/obj.H;
            D=T;

            Lp=L*obj.AeroM;
            Dp=D*obj.AeroM;


            %Center of mass
            coordCOM=zeros(obj.n_d,1);
            sumX=0;
            sumY=0;
            sumZ=0;

            % Pes(3)=Pes(3)-Wm/(2);
            % Pes(6)=Pes(6)-Wm/(2);

            massa=-Pes/9.81;

            massa(3)=massa(3)+obj.Wm/(2*9.81);
            massa(6)=massa(6)+obj.Wm/(2*9.81);

            for i=1:obj.n
                sumX=sumX+obj.x(i,1)*massa(3*i);
                sumY=sumY+obj.x(i,2)*massa(3*i);
                sumZ=sumZ+obj.x(i,3)*massa(3*i);
            end




            coordCOM(1)=sumX/(mTotal+obj.Wm/9.81);
            coordCOM(2)=sumY/(mTotal+obj.Wm/9.81);
            coordCOM(3)=sumZ/(mTotal+obj.Wm/9.81);

            Tp=(3/5*Lp*coordCOM(1)+Lp/5*(coordCOM(1)-obj.W)-Lp/5*(2*obj.W-coordCOM(1))-Dp*(obj.H-coordCOM(3)))/coordCOM(3);

            aX=(Tp-Dp)/(mTotal+obj.Wm/9.81);
            aZ=(Lp-obj.Wm-Ptotal)/(mTotal+obj.Wm/9.81);



            Fin=zeros(obj.n_dof,1);

            if obj.AeroM==1
                Fin=zeros(obj.n_dof,1);

                Fdata = [1 1 T/2 ; 2 4 T/2 ; 1 3 -obj.Wm/2 ; 2 6 -obj.Wm/2 ; 7 21 L/5 ; 3 9 L/5 ; 4 12 L/5 ; 5 15 L/5 ; 6 18 L/5 ;
                    7 19 -D/5 ; 3 7 -D/5 ; 4 10 -D/5 ; 5 13 -D/5 ; 6 16 -D/5 ]; %en globals
            else
                for i=1:obj.n
                    Fin(3*i,:)=-massa(3*i)*aZ;
                end


                Fdata = [1 1 Tp/2 ; 2 4 Tp/2 ; 1 3 -obj.Wm/2 ; 2 6 -obj.Wm/2 ; 7 21 Lp/5 ; 3 9 Lp/5 ; 4 12 Lp/5 ; 5 15 Lp/5 ; 6 18 Lp/5 ;
                    7 19 -Dp/5 ; 3 7 -Dp/5 ; 4 10 -Dp/5 ; 5 13 -Dp/5 ; 6 16 -Dp/5 ];
            end

            c=size(Fdata,1);

            for i=1:c
                F_ext(Fdata(i,2),1) = Fdata(i,3);
            end

            Fext=F_ext+Pes+Fin;


        end
    end



    methods (Access = private)

        function init(obj, cParams)
            obj.n = cParams.dimensionalData.n;
            obj.n_d = cParams.dimensionalData.n_d;
            obj.n_el = cParams.dimensionalData.n_el;
            obj.n_dof = cParams.dimensionalData.n_dof;
            obj.x = cParams.preprocessData.coor;
            obj.Tn = cParams.preprocessData.nodalConnec;
            obj.mat = cParams.preprocessData.matProp;
            obj.Tmat = cParams.preprocessData.matConnec;
            obj.Wm = 9.81*cParams.hangingMass;
            obj.AeroM = cParams.aeroMultiplier;

            obj.W = cParams.geometricalData.W;
            obj.H = cParams.geometricalData.H;
            obj.D1 = cParams.geometricalData.D1;
            obj.d1 = cParams.geometricalData.d1;
            obj.D2 = cParams.geometricalData.D2;
        end

    end


end

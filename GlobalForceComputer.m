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
        hWeight
        AeroM        
        W
        H
        D1
        d1
        D2
    end
    
    properties (Access = private)
        wVec
        Ptotal
        coordCOM
        aeroForces
        aeroMultForces
        aZ
        Fext
    end
    
    methods (Access = public)
        
        function obj = GlobalForceComputer(cParams)
            obj.init(cParams)
        end
        
        function externalForces = compute(obj)            
            obj.computeNodalWeight();            
            obj.computeCenterOfMass();            
            obj.computeAeroForces();            
            obj.computeAccelerations();            
            obj.computeForces();
            externalForces = obj.Fext;
        end
    end
        
    methods (Access = private)
        
        function init(obj, cParams)
            obj.n       = cParams.n;
            obj.n_d     = cParams.n_d;
            obj.n_el    = cParams.n_el;
            obj.n_dof   = cParams.n_dof;
            obj.x       = cParams.mesh.coor;
            obj.Tn      = cParams.mesh.nodalConnec;
            obj.mat     = cParams.materialData.matProp;
            obj.Tmat    = cParams.materialData.matConnec;
            obj.hWeight = cParams.hMass*9.81;
            obj.AeroM   = cParams.AeroM;            
            obj.W       = cParams.W;
            obj.H       = cParams.H;
            obj.D1      = cParams.D1;
            obj.d1      = cParams.d1;
            obj.D2      = cParams.D2;
        end
                
        function computeNodalWeight(obj)
            obj.Ptotal = 0;
            Pes = zeros(obj.n_dof,1);            
            for iElem = 1:obj.n_el
                nodeA   = obj.Tn(iElem,1);
                nodeB   = obj.Tn(iElem,2);                
                l_e     = obj.computeElementLength(nodeA,nodeB);
                V       = obj.computeElementVolume(iElem, l_e);                
                obj.Ptotal = obj.Ptotal + V*obj.mat(obj.Tmat(iElem),3)*9.81;
                Pes(3*nodeA) = Pes(3*nodeA) - V*obj.mat(obj.Tmat(iElem),3)*9.81/2;
                Pes(3*nodeB) = Pes(3*nodeB) - V*obj.mat(obj.Tmat(iElem),3)*9.81/2;
            end                      
            obj.wVec = Pes;
        end
        
        function elementLength = computeElementLength(obj,nodeA,nodeB)            
            s.x1 = obj.x(nodeA,1);
            s.y1 = obj.x(nodeA,2);
            s.z1 = obj.x(nodeA,3);
            s.x2 = obj.x(nodeB,1);
            s.y2 = obj.x(nodeB,2);
            s.z2 = obj.x(nodeB,3);
            elementLength = sqrt((s.x2-s.x1)^2 + (s.y2-s.y1)^2+(s.z2-s.z1)^2);
        end
        
        function elementVolume = computeElementVolume(obj, e, elementLength)
            if obj.Tmat(e) == 1
                elementVolume=pi*(obj.D1/2)^2*elementLength-pi*(obj.d1/2)^2*elementLength;
            else
                elementVolume=(obj.D2/2)^2*pi*elementLength;
            end
        end
        
        function computeAeroForces(obj)
            M_PES_7 = 0;
            for i = 1:obj.n
                M_PES_7 = M_PES_7 + obj.wVec(3*i)*(obj.x(7,1) - obj.x(i,1));
            end            
            aeroF.L = obj.hWeight+obj.Ptotal;
            aeroF.T = (obj.hWeight*obj.W+2/5*aeroF.L*obj.W+M_PES_7)/obj.H;
            aeroF.D = aeroF.T;            
            aeroFp.L = aeroF.L*obj.AeroM;
            aeroFp.D = aeroF.D*obj.AeroM;            
            aeroFp.T = (3/5*aeroFp.L*obj.coordCOM(1)+aeroFp.L/5*(obj.coordCOM(1)-obj.W)-aeroFp.L/5*(2*obj.W-obj.coordCOM(1))-aeroFp.D*(obj.H-obj.coordCOM(3)))/obj.coordCOM(3);
            obj.aeroForces = aeroF;
            obj.aeroMultForces = aeroFp;
        end
        
        function computeCenterOfMass(obj)
            mTotal  = obj.Ptotal/9.81;
            COM     = zeros(obj.n_d,1);
            sumX    = 0;
            sumY    = 0;
            sumZ    = 0;            
            massVec     = -obj.wVec/9.81;            
            massVec(3)  = massVec(3) + obj.hWeight/(2*9.81);
            massVec(6)  = massVec(6) + obj.hWeight/(2*9.81);            
            for i=1:obj.n
                sumX    = sumX+obj.x(i,1)*massVec(3*i);
                sumY    = sumY+obj.x(i,2)*massVec(3*i);
                sumZ    = sumZ+obj.x(i,3)*massVec(3*i);
            end            
            COM(1)  = sumX/(mTotal+obj.hWeight/9.81);
            COM(2)  = sumY/(mTotal+obj.hWeight/9.81);
            COM(3)  = sumZ/(mTotal+obj.hWeight/9.81);
            obj.coordCOM = COM;
        end
        
        function computeAccelerations(obj)
            mTotal = obj.Ptotal;
            obj.aZ = (obj.aeroMultForces.L-obj.hWeight-obj.Ptotal)/(mTotal+obj.hWeight/9.81);
        end
        
        function computeForces(obj)
            F_ext   = zeros(obj.n_dof,1);
            Fin     = zeros(obj.n_dof,1);
            massVec = obj.wVec/9.81;            
            if obj.AeroM == 1
                aero = obj.aeroForces;
            else
                for i = 1:obj.n
                    Fin(3*i,:) = -massVec(3*i)*obj.aZ;
                end
                aero = obj.aeroMultForces;
            end            
            Fdata = [   1 1 aero.T/2 ;
                        2 4 aero.T/2 ;
                        1 3 -obj.hWeight/2 ;
                        2 6 -obj.hWeight/2 ;
                        7 21 aero.L/5 ;
                        3 9 aero.L/5 ;
                        4 12 aero.L/5 ;
                        5 15 aero.L/5 ;
                        6 18 aero.L/5 ;
                        7 19 -aero.D/5 ;
                        3 7 -aero.D/5 ;
                        4 10 -aero.D/5 ;
                        5 13 -aero.D/5 ;
                        6 16 -aero.D/5];            
            c = size(Fdata,1);            
            for i = 1:c
                F_ext(Fdata(i,2),1) = Fdata(i,3);
            end
            obj.Fext = F_ext+obj.wVec+Fin;
        end        
    end  
end

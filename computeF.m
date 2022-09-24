function [Fext,aX,aZ] = computeF(n_d,n_el,n_dof,W,H,D1,d1,D2,x,Tmat,mat,Tn,n,Wm,AeroM)
%--------------------------------------------------------------------------
% The function takes as inputs:
%   - Dimensions:  n_i         Number of DOFs per node
%                  n_dof       Total number of DOFs
%   - Fdata  External nodal forces [Nforces x 3]
%            Fdata(k,1) - Node at which the force is applied
%            Fdata(k,2) - DOF (direction) at which the force acts
%            Fdata(k,3) - Force magnitude in the corresponding DOF
%--------------------------------------------------------------------------
% It must provide as output:
%   - Fext  Global force vector [n_dof x 1]
%            Fext(I) - Total external force acting on DOF I
%--------------------------------------------------------------------------
% Hint: Use the relation between the DOFs numbering and nodal numbering to
% determine at which DOF in the global system each force is applied.
F_ext=zeros(n_dof,1);


Pes=zeros(n_dof,1);
massa=zeros(n_dof,1);

x1=zeros(n_el,1);
y1=zeros(n_el,2);
z1=zeros(n_el,3);

x2=zeros(n_el,1);
y2=zeros(n_el,2);
z2=zeros(n_el,3);
Ptotal=0;

for e=1:n_el
    x1(e) = x(Tn(e,1),1);
    y1(e) = x(Tn(e,1),2);
    z1(e) = x(Tn(e,1),3);
    x2(e) = x(Tn(e,2),1);
    y2(e) = x(Tn(e,2),2);
    z2(e) = x(Tn(e,2),3);
  
    
    l(e)=sqrt((x2(e)-x1(e))^2 + (y2(e)-y1(e))^2+(z2(e)-z1(e))^2);
    n1=Tn(e,1);
    n2=Tn(e,2);
    
    if Tmat(e) == 1
        
        V=pi*(D1/2)^2*l(e)-pi*(d1/2)^2*l(e);    
        
    else
            
        V=(D2/2)^2*pi*l(e);
        
    end
        
    Ptotal=Ptotal+V*mat(Tmat(e),3)*9.81;
    
    Pes(3*n1)=Pes(3*n1)-V*mat(Tmat(e),3)*9.81/2;
    Pes(3*n2)=Pes(3*n2)-V*mat(Tmat(e),3)*9.81/2;
    
end
% 
% massa=-Pes/9.81;
mTotal=Ptotal/9.81;


M_PES_7=0;

for i=1:n
    M_PES_7=M_PES_7+Pes(3*i)*[x(7,1)-x(i,1)];
end


L=Wm+Ptotal;
T=(Wm*W+2/5*L*W+M_PES_7)/H;
D=T;

Lp=L*AeroM;
Dp=D*AeroM;


%Center of mass
coordCOM=zeros(n_d,1);
sumX=0;
sumY=0;
sumZ=0;

% Pes(3)=Pes(3)-Wm/(2);
% Pes(6)=Pes(6)-Wm/(2);

massa=-Pes/9.81;

massa(3)=massa(3)+Wm/(2*9.81);
massa(6)=massa(6)+Wm/(2*9.81);

for i=1:n
    sumX=sumX+x(i,1)*massa(3*i);
    sumY=sumY+x(i,2)*massa(3*i);
    sumZ=sumZ+x(i,3)*massa(3*i);
end




coordCOM(1)=sumX/(mTotal+Wm/9.81);
coordCOM(2)=sumY/(mTotal+Wm/9.81);
coordCOM(3)=sumZ/(mTotal+Wm/9.81);

Tp=(3/5*Lp*coordCOM(1)+Lp/5*(coordCOM(1)-W)-Lp/5*(2*W-coordCOM(1))-Dp*(H-coordCOM(3)))/coordCOM(3);

aX=(Tp-Dp)/(mTotal+Wm/9.81);
aZ=(Lp-Wm-Ptotal)/(mTotal+Wm/9.81);




if AeroM==1
    Fin=zeros(n_dof,1);
    
    Fdata = [1 1 T/2 ; 2 4 T/2 ; 1 3 -Wm/2 ; 2 6 -Wm/2 ; 7 21 L/5 ; 3 9 L/5 ; 4 12 L/5 ; 5 15 L/5 ; 6 18 L/5 ;
        7 19 -D/5 ; 3 7 -D/5 ; 4 10 -D/5 ; 5 13 -D/5 ; 6 16 -D/5 ]; %en globals
else
    for i=1:n
        Fin(3*i,:)=-massa(3*i)*aZ;         
    end
    
        
    Fdata = [1 1 Tp/2 ; 2 4 Tp/2 ; 1 3 -Wm/2 ; 2 6 -Wm/2 ; 7 21 Lp/5 ; 3 9 Lp/5 ; 4 12 Lp/5 ; 5 15 Lp/5 ; 6 18 Lp/5 ;
        7 19 -Dp/5 ; 3 7 -Dp/5 ; 4 10 -Dp/5 ; 5 13 -Dp/5 ; 6 16 -Dp/5 ]; 
end

c=size(Fdata,1);

for i=1:c
    F_ext(Fdata(i,2),1) = Fdata(i,3);
end
    
Fext=F_ext+Pes+Fin;


end

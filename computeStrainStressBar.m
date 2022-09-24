function [eps,sig] = computeStrainStressBar(n_d,n_i,n_nod,n_el,u,Td,x,Tn,mat,Tmat)
%--------------------------------------------------------------------------
% The function takes as inputs:
%   - Dimensions:  n_d        Problem's dimensions
%                  n_el       Total number of elements
%   - u     Global displacement vector [n_dof x 1]
%            u(I) - Total displacement on global DOF I
%   - Td    DOFs connectivities table [n_el x n_el_dof]
%            Td(e,i) - DOF i associated to element e
%   - x     Nodal coordinates matrix [n x n_d]
%            x(a,i) - Coordinates of node a in the i dimension
%   - Tn    Nodal connectivities table [n_el x n_nod]
%            Tn(e,a) - Nodal number associated to node a of element e
%   - mat   Material properties table [Nmat x NpropertiesXmat]
%            mat(m,1) - Young modulus of material m
%            mat(m,2) - Section area of material m
%   - Tmat  Material connectivities table [n_el]
%            Tmat(e) - Material index of element e
%--------------------------------------------------------------------------
% It must provide as output:
%   - eps   Strain vector [n_el x 1]
%            eps(e) - Strain of bar e
%   - sig   Stress vector [n_el x 1]
%            sig(e) - Stress of bar e
%--------------------------------------------------------------------------

c=size(Td,2);
n_el_dof=n_d*n_nod;

x1=zeros(n_el,1);
y1=zeros(n_el,2);
z1=zeros(n_el,3);

x2=zeros(n_el,1);
y2=zeros(n_el,2);
z2=zeros(n_el,3);

l=zeros(n_el,1);

Re=zeros(n_nod,n_el_dof,n_el);
u_e=zeros(c,n_el);
u_p=zeros(n_nod,n_el);
eps=zeros(n_el,1);
sig=zeros(n_el,1);



for e=1:n_el
        E=mat(Tmat(e),1);
    A=mat(Tmat(e),2);
    x1(e) = x(Tn(e,1),1);
    y1(e) = x(Tn(e,1),2);
    z1(e) = x(Tn(e,1),3);
    x2(e) = x(Tn(e,2),1);
    y2(e) = x(Tn(e,2),2);
    z2(e) = x(Tn(e,2),3);
  
    
    l(e)=sqrt((x2(e)-x1(e))^2 + (y2(e)-y1(e))^2+(z2(e)-z1(e))^2);
    
    Re(:,:,e)=1/l(e)*[x2(e)-x1(e) y2(e)-y1(e) z2(e)-z1(e) 0 0 0;
                        0 0 0 x2(e)-x1(e) y2(e)-y1(e) z2(e)-z1(e)];

    
      

    
    for i=1:n_i*n_nod
        
        I=Td(e,i);
        u_e(i,e)=u(I);
        
    end
    
         u_p(:,e)=Re(:,:,e)*u_e(:,e);

         eps(e,1)=(1/l(e))*[-1 1]*u_p(:,e);

         sig(e,1)=E*eps(e,1);
    

end


end
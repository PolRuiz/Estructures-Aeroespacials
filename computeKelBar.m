function Kel = computeKelBar(n_d,n_el,x,Tn,mat,Tmat)
%--------------------------------------------------------------------------
% The function takes as inputs:
%   - Dimensions:  n_d        Problem's dimensions
%                  n_el       Total number of elements
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
%   - Kel   Elemental stiffness matrices [n_el_dof x n_el_dof x n_el]
%            Kel(i,j,e) - Term in (i,j) position of stiffness matrix for element e
%--------------------------------------------------------------------------

n_nod=size(Tn,2);
n_el_dof=n_d*n_nod;

x1=zeros(n_el,1);
y1=zeros(n_el,2);
z1=zeros(n_el,3);

x2=zeros(n_el,1);
y2=zeros(n_el,2);
z2=zeros(n_el,3);


Re=zeros(n_nod,n_el_dof);
Kep=zeros(n_nod,n_nod);
Ke=zeros(n_el_dof, n_el_dof);

Kel=zeros(n_el_dof, n_el_dof,n_el);

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
    
    Re=1/l(e)*[x2(e)-x1(e) y2(e)-y1(e) z2(e)-z1(e) 0 0 0;
                        0 0 0 x2(e)-x1(e) y2(e)-y1(e) z2(e)-z1(e)];
    
    
    Kep=(A*E)/l(e)*[1 -1;-1 1];
    
    Ke=transpose(Re)*Kep*Re;
    
    Kel(:,:,e)=Ke(:,:);
end


end
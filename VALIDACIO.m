%-------------------------------------------------------------------------%
% ASSIGNMENT 02
%-------------------------------------------------------------------------%
% Date:
% Author/s:
%

clear;
close all;

%% INPUT DATA
H = 0.9;
W = 0.85;
B = 3.2;
D1 = 18/1000; %m
d1 = 7.5/1000; %m
D2 = 3/1000; %m

% Mass
M = 150;

% Other
g = 9.81;

%Forces
F=100;

%Multiplicador F.Aerodinamiques

AeroM = 1;
% Geometric data

L=2;



% Other
g = 9.81;

%% PREPROCESS

% Nodal coordinates matrix 
%  x(a,j) = coordinate of node a in the dimension j
x = [%     X      Y      Z
         0,  0,     0; % (1)
         0,   L,     0; % (2)
         0,     L,     L; % (3)
         0,     0,     L; % (4)  
         L,     L/2,     L/2; % (5)
];

% Nodal connectivities  
%  Tnod(e,a) = global nodal number associated to node a of element e
Tn = [1 5; 2 5; 3 5; 4 5];

% Material properties matrix
%  mat(m,1) = Young modulus of material m
%  mat(m,2) = Section area of material m
%  mat(m,3) = Density of material m
%  --more columns can be added for additional material properties--
mat = [% Young M.        Section A.    Density   Inertia             
                147000000000,                pi*(D2/2)^2,           950,      pi/32*D2^4;% Material (2)
];

% Material connectivities
%  Tmat(e) = Row in mat corresponding to the material associated to element e 
Tmat = [1 1 1 1];
 


% Fix nodes matrix creation
%  fixNod(k,1) = node at which some DOF is prescribed
%  fixNod(k,2) = DOF prescribed
%  fixNod(k,3) = prescribed displacement in the corresponding DOF (0 for fixed)

fixNod = [1 1 0; 1 2 0; 1 3 0; 2 4 0; 2 5 0; 2 6 0; 3 7 0; 3 8 0; 3 9 0; 4 10 0; 4 11 0; 4 12 0];

%% SOLVER

% Dimensions
n_d = size(x,2);              % Number of dimensions
n_i = n_d;                    % Number of DOFs for each node
n = size(x,1);                % Total number of nodes
n_dof = n_i*n;                % Total number of degrees of freedom
n_el = size(Tn,1);            % Total number of elements
n_nod = size(Tn,2);           % Number of nodes for each element
n_el_dof = n_i*n_nod;         % Number of DOFs for each element 

% Computation of the DOFs connectivities
Td = connectDOFs(n_el,n_nod,n_i,Tn);

% Computation of element stiffness matrices
Kel = computeKelBar(n_d,n_el,x,Tn,mat,Tmat);

% Global matrix assembly
KG = assemblyKG(n_el,n_el_dof,n_dof,Td,Kel);

% Global force vector assembly
Fext = computeF2(n_dof,F);

% Apply conditions 
[vL,vR,uR] = applyCond(n_i,n_dof,fixNod);

% System resolution
[u,R] = solveSys(vL,vR,uR,KG,Fext);

% Compute strain and stresses
[eps,sig] = computeStrainStressBar(n_d,n_i,n_nod,n_el,u,Td,x,Tn,mat,Tmat);


% Buckling

[sig_cr] = buckling(n_d,n_el,x,mat,Tn,Tmat,sig);


%% POSTPROCESS

% Plot deformed structure with stress of each bar
scale = 10; % Adjust this parameter for properly visualizing the deformation
plotBarStress3D(x,Tn,u,sig,scale);
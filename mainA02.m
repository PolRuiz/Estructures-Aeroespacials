%-------------------------------------------------------------------------%
% ASSIGNMENT 2
%-------------------------------------------------------------------------%
% Date: 09/03/2022
% Author/s: Albert Servitje Roca i Pol Ruiz Celada
%
clc;
clear;
close all;

%% INPUT DATA

% Geometric data
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
Wm=M*g;

%Multiplicador F.Aerodinamiques

AeroM = 1;



%% PREPROCESS

% Nodal coordinates matrix 
%  x(a,j) = coordinate of node a in the dimension j
x = [%     X      Y      Z
         2*W,  -W/2,     0; % (1)
         2*W,   W/2,     0; % (2)
         2*W,     0,     H; % (3)
           0,     0,     H; % (4)
           0,    -B,     H; % (5)
           0,     B,     H; % (6)
           W,     0,     H; % (7)
];

% Nodal connectivities  
%  Tnod(e,a) = global nodal number associated to node a of element e
Tn = [1 2;2 3;1 3;3 5;3 6;4 6;6 7;4 5;5 7;3 7;4 7;1 5;2 6;1 7;2 7;1 4;2 4];

% Material properties matrix
%  mat(m,1) = Young modulus of material m
%  mat(m,2) = Section area of material m
%  mat(m,3) = Density of material m
%  --more columns can be added for additional material properties--
mat = [% Young M.        Section A.    Density   Inertia
                75000000000,                pi*(D1/2)^2-pi*(d1/2)^2,           3350,  pi/32*(D1^4-d1^4);% Material (1)
                147000000000,                pi*(D2/2)^2,           950,      pi/32*D2^4;% Material (2)
];

% Material connectivities
%  Tmat(e) = Row in mat corresponding to the material associated to element e 
Tmat = [1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2];
 

% External force matrix creation
%  Fdata(k,1) = node at which the force is applied
%  Fdata(k,2) = DOF (direction) at which the force is applied
%  Fdata(k,3) = force magnitude in the corresponding DOF
% Fdata = [1 1 T/2 ; 2 4 T/2 ; 1 3 -W/2 ; 2 6 -W/2 ; 7 21 L/5 ; 3 9 L/5 ; 4 12 L/5 ; 5 15 L/5 ; 6 18 L/5 ;
%     7 19 -D/5 ; 3 7 -D/5 ; 4 10 -D/5 ; 5 13 -D/5 ; 6 16 -D/5 ]; %en globals

% Fix nodes matrix creation
%  fixNod(k,1) = node at which some DOF is prescribed
%  fixNod(k,2) = DOF prescribed
%  fixNod(k,3) = prescribed displacement in the corresponding DOF (0 for fixed)
fixNod = [4 10 0; 4 11 0; 4 12 0; 5 13 0; 5 15 0; 7 21 0];


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
% KG = assemblyKG(n_el,n_el_dof,n_dof,Td,Kel);
%%
stiffnessComputer = GlobalStiffnessMatrixComputer(Kel,Tn,n,n_i);
KG = stiffnessComputer.obtainStiffnessMatrix();
%%


% Global force vector assembly
[Fext,aX,aZ] = computeF(n_d,n_el,n_dof,W,H,D1,d1,D2,x,Tmat,mat,Tn,n,Wm,AeroM);

% Apply conditions 
[vL,vR,uR] = applyCond(n_i,n_dof,fixNod);

% System resolution
solverType = 'Iterative'; %'Iterative' or 'Direct'
[uL,u,R] = solveSys(vL,vR,uR,KG,Fext,solverType);

% Compute strain and stresses
[eps,sig] = computeStrainStressBar(n_d,n_i,n_nod,n_el,u,Td,x,Tn,mat,Tmat);


% Buckling
[sig_cr] = buckling(n_d,n_el,x,mat,Tn,Tmat,sig);


%% POSTPROCESS

% Plot deformed structure with stress of each bar
scale = 30; % Adjust this parameter for properly visualizing the deformation
plotBarStress3D(x,Tn,u,sig,scale);

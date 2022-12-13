clc;
clear;
close all;



%% INPUT DATA

% Mass
M = 150;

% Other
g = 9.81;

%Forces
Wm=M*g;

%Multiplicador F.Aerodinamiques
AeroM = 1;


%%
problem = DataCreator();
Data = problem.createData();

x = Data.preprocessData.coor;
Tn = Data.preprocessData.nodalConnec;
mat = Data.preprocessData.matProp;
Tmat = Data.preprocessData.matConnec;
fixNod = Data.preprocessData.fixNode;

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

%%
stiffnessComputer = GlobalStiffnessMatrixComputer(Kel,Tn,n,n_i);
KG = stiffnessComputer.obtainStiffnessMatrix();
%%


% Global force vector assembly
[Fext,aX,aZ] = computeF(n_d,n_el,n_dof,problem,n,Wm,AeroM);

% Apply conditions 
[vL,vR,uR] = applyCond(n_i,n_dof,fixNod);

% System resolution
solverType = 'Direct'; %'Iterative' or 'Direct'
[uL,u,R] = solveSys(vL,vR,uR,KG,Fext,solverType);

% Compute strain and stresses
[eps,sig] = computeStrainStressBar(n_d,n_i,n_nod,n_el,u,Td,x,Tn,mat,Tmat);


% Buckling
[sig_cr] = buckling(n_d,n_el,x,mat,Tn,Tmat,sig);


%% POSTPROCESS

% Plot deformed structure with stress of each bar
scale = 30; % Adjust this parameter for properly visualizing the deformation
plotBarStress3D(x,Tn,u,sig,scale);

%% TESTS

runtests('tests.m')

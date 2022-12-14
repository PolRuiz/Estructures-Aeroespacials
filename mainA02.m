clc;
clear;
close all;

%%
s.Mass = 150; %kg
s.AeroM = 1;

hangGlider = TrussStructure(s);

hangGlider.createData();

hangGlider.computeDOFConnectivity();

hangGlider.computeElementStiffnes();

hangGlider.computeStiffnessMatrix();
KG = hangGlider.KG;

hangGlider.computeGlobalForceVector();
Fext = hangGlider.Fext;

hangGlider.applyBoundaryCond();
vL = hangGlider.vectorialData.vL;
vR = hangGlider.vectorialData.vR;


% System resolution

solverType = 'Direct'; %'Iterative' or 'Direct'
[uL,u,R] = solveSys(hangGlider, KG, Fext, solverType);

hangGlider.computeStrainStress(u, uL);


%% POSTPROCESS

% Plot deformed structure with stress of each bar
scale = 30; % Adjust this parameter for properly visualizing the deformation
plotBarStress3D(hangGlider,scale);

%% TESTS

runtests('tests.m')

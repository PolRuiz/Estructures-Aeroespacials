clc;clear;close all;

dataCreator = DataCreator();
problemData = dataCreator.createData();
hangGlider  = TrussStructure(problemData);
hangGlider.computeProblem();
runtests('tests.m')

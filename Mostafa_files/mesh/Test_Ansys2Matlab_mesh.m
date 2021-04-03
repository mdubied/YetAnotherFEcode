clear all
close all
meshfile='MeshFreqDividerFull_QUAD4.dat';

meshtype='QUAD4';

model=Ansys2Matlab_mesh(meshfile,meshtype);

PlotMesh(model.nodes(:,2:3),model.elements(:,2:end),0)



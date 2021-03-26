clear all
close all
meshfile='MeshBeam2D_QUAD8_metric.dat';
meshtype='QUAD8';

model=Ansys2Matlab_mesh(meshfile,meshtype);

PlotMesh(model.nodes(:,2:3),model.elements(:,2:end),0)



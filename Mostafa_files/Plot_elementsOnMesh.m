
figure
hold on
PlotMesh(nodes,elements(:,1:4),0)

elementHyper=[];
elementWeightfnnls=[];
for ii=1:size(x_fnnls,1)
    if x_fnnls(ii)~=0
        elementHyper=[elementHyper ii];
        elementWeightfnnls=[elementWeightfnnls x_fnnls(ii)];
    end
end


PlotMesh(nodes,elements([elementHyper],1:4),0)

elementHypernnls=[];
elementWeightnnls=[];
for ii=1:size(x_nnls,1)
    if x_nnls(ii)~=0
        elementHypernnls=[elementHypernnls ii];
        elementWeightnnls=[elementWeightnnls x_nnls(ii)];
    end
end

PlotMesh(nodes,elements([elementHypernnls],1:4),0)

elementHypernnlsJ=[];
elementWeightnnlsJ=[];
for ii=1:size(xi_jain_d,1)
    if xi_jain_d(ii)~=0
        elementHypernnlsJ=[elementHypernnlsJ ii];
        elementWeightnnlsJ=[elementWeightnnlsJ xi_jain_d(ii)];
    end
end

PlotMesh(nodes,elements([elementHypernnlsJ],1:4),0)
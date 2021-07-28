[UU SS VV]=svd(  full(bForce)) ;

y1SS=diag(SS);
xSS=linspace(1,size(y1SS,1),size(y1SS,1));
figure 
hold on
plot(xSS,normalize(y1SS,'range'),'+-')

for mod=1:5
    % mod = 4;
    elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
    figure('units','normalized','position',[.2 .1 .6 .8])
    
    PlotMesh(nodes, elementPlot, 0);
    v1 = reshape(UU(:,mod), 2, []).';
    PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor',1e-4);
   % title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0(mod),3) ' Hz'])
    
end
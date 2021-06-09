W=[];
W_names=[];
eta=VMs'*TI_lin.Solution.u;
for ii=1:size(eta,1)
    for jj=ii:size(eta,1)
       W_names=[ii ;jj];
       MMI= max(eta(ii,:).*eta(jj,:));
       W=[W [W_names ; MMI]];
    end
end

figure()

hold on

scatterbar3(W(1,:),W(2,:),W(3,:),1)
scatter3(W(1,:),W(2,:),W(3,:),'*')
grid on ;box on;
xticks(unique(W(1,:)))
yticks(unique(W(2,:)))

sortedW=sortrows(W',3,'descend')';

sortedMDs=[];
for ii=1:size(sortedW,2)
    I_index=sortedW(1,ii);
    J_index=sortedW(2,ii);
    MD_ii=find(and(MDs_names(:,1)== I_index,MDs_names(:,2)== J_index));
    sortedMDs=[sortedMDs MD_ii];
end

% MD_used=[];
% 
% MD1=find(and(MDs_names(:,1)== 1,MDs_names(:,2)== 3));
% MD2=find(and(MDs_names(:,1)== 1,MDs_names(:,2)== 1));
% MD3=find(and(MDs_names(:,1)== 3,MDs_names(:,2)== 3));
% MD4=find(and(MDs_names(:,1)== 1,MDs_names(:,2)== 6));


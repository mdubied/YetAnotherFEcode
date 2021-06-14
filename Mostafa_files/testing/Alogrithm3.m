%% Algorithm 3 

n_V=3;
%model superposition
V_H=VMs(:,1:n_V)
for ii = 1 : size(V_H, 2)
             V_H(:,ii) = V_H(:,ii) / (V_H(:,ii)'*(BeamAssembly.DATA.M)*V_H(:,ii));
end

Lin_sol=TI_lin.Solution.u;
Lin_sol_snap=[];
     for ii=1:size(Lin_sol,2)
         if rem(ii,1)==0
         Lin_sol_snap=[Lin_sol_snap Lin_sol(:,ii)];
         end
     end
eta=V_H'*Lin_sol_snap;
%MD calculations
[Th MDnames] = modal_derivatives(BeamAssembly, elements, V_H);

Theta = QM_Theta_from_SMDs(Th, MDnames);

%perform MD selection
V_ECSW=[V_H Th];
%QM uplifting
uu = zeros(size(V_H,1), size(eta, 2));
for tt = 1 : size(eta, 2)
    uu(:,tt) = einsum('Iij,iJ,jK->IJK', Theta, eta(:,tt), eta(:,tt));
end
u_lin_ECSW = V_H * eta + 1/2*uu;

%Construct Gb
qq=(V_ECSW'*V_ECSW)^-1*V_ECSW'*u_lin_ECSW;

% nt=size(qq,2);
% G=zeros(size(V_ECSW,2)*nt,BeamAssembly.Mesh.nElements);
% b=zeros(size(V_ECSW,2)*nt,1)
%  for ii=1:nt
%      [~,Fv]=BeamAssembly.tangent_stiffness_and_force(V_ECSW*qq(:,ii));
%      G(ii,:)=Fv';
%      b(ii)=sum(Fv);
%  end
tic
[G,b]=BeamReducedAssembly.constructGb(qq);
toc

%fNNLS
tic
[xii,ww]=fnnls(G'*G,G'*b,0.001);
toc

disp('nnls')
tic;[xx1,ww1,info]=nnls(full(G),full(b),struct('Accy',2));toc
disp(['Min x  ', num2str(min(xx1)), '    Max w  ', num2str(max(ww1(xx1==0)))])
disp(' ')
disp('lsqnoneg ')
options = optimset('TolX',0.01);
tic;x_lsg=lsqnonneg(G,b);toc
disp(['Min x  ', num2str(min(x_lsg)), '    Difference x-x1  ', num2str(max(abs(xx1-x_lsg)))])
disp(' ')


function Theta = QM_Theta_from_SMDs(SMDs, names)

    n = size(SMDs,1);
    m = max(names(:));
    I = names(:,1);
    J = names(:,2);

    Theta = zeros(n,m,m);
    for k = 1:size(SMDs,2)
        Theta(:,I(k),J(k)) = SMDs(:,k);
        Theta(:,J(k),I(k)) = SMDs(:,k);
    end

end

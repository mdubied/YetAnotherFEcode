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
         if rem(ii,10)==0
         Lin_sol_snap=[Lin_sol_snap Lin_sol(:,ii)];
         end
     end
eta=V_H'*Lin_sol_snap;
%MD calculations
[Th MDnames] = modal_derivatives(BeamAssembly, elements, V_H);
% for ii = 1:size(Th,2)
%     Th(:,ii) = Th(:,ii)/max(sqrt(sum(Th(:,ii).^2,2)));
%  end
Theta = QM_Theta_from_SMDs(Th, MDnames);

 [DSs, DSnames] = defect_sensitivities(BeamAssembly, elements, V_H, U, ...
     FORMULATION);
%  for ii = 1:size(DS,2)
%     DS(:,ii) = DS(:,ii)/max(sqrt(sum(DS(:,ii).^2,2)));
%  end
Xi = DS_Xi_QM(DSs, DSnames);

%perform MD selection
V_ECSW=[VMs(:,1:n_V) Th DSs];

for ii = 1 : size(V_ECSW, 2)
             V_ECSW(:,ii) = V_ECSW(:,ii) / (V_ECSW(:,ii)'*(BeamAssembly.DATA.M)*V_ECSW(:,ii));
end
%QM uplifting
uu = zeros(size(V_H,1), size(eta, 2));

for tt = 1 : size(eta, 2)
    uu(:,tt) = einsum('Iij,iJ,jK->IJK', Theta, eta(:,tt), eta(:,tt));
end

uud = zeros(size(V_H,1), size(eta, 2));

for tt = 1 : size(eta, 2)
    for ii= 1:size(xi,2)
    uud(:,tt) =uud(:,tt)+ einsum('Iij,iJ,jK->IJK', Xi, eta(:,tt), xi(:,ii));
end
end
u_lin_ECSW = V_H * eta + 1/2*uu+uud;

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
[G,b]=BeamDefectedReducedAssembly.constructGb(qq);
toc

%fNNLS
tic
[x_fnnls,w_fnnls]=fnnls(G'*G,G'*b,1.1722e+06);
toc

disp('nnls')
tic;[x_nnls,w_nnls,info]=nnls(full(G),full(b),struct('Accy',1,'Tol',10e-10));toc
disp(['Min x  ', num2str(min(x_nnls)), '    Max w  ', num2str(max(w_nnls(x_nnls==0)))])
disp(' ')
disp('lsqnoneg ')
 options = optimset('TolX',1.1722e+06);
tic;x_lsq=lsqnonneg(G,b,options);toc
% disp(['Min x  ', num2str(min(x_lsg)), '    Difference x-x1  ', num2str(max(abs(xx1-x_lsg)))])
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
function Xi = DS_Xi_QM(DSs, names)

    n = size(DSs,1);
    m = max(names(:,1));
    md = max(names(:,2));

    I = names(:,1);
    J = names(:,2);

    Xi = zeros(n,m,md);
    for k = 1:size(DSs,2)
        Xi(:,I(k),J(k)) = DSs(:,k);
    end

end
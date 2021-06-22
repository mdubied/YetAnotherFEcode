%to be deleted
numofELEMENTS=[];
figure 
hold on
for kkk=[5]
LROMHd=[];
LtimeHd=[];
LdurationHd=[];
NLROMHd=[];
NLtimeHd=[];
NLdurationHd=[];
f0c_ROMHd={};

nfsense = find_node(Lx/2,Ly/2,[],nodes); % node where to see the results
dof = get_index(nfsense, myMesh.nDOFPerNode );
dof=dof(2);

countHd=0;
VbModeHd=[];
MDeriHd=[];
DSentiHd=[];
for m=3;% m=[2 3 7 10]
    for n=m*(m+1)/2
     for k=m
        
        t=m+n+k;
      
%%
        
        if n==0
            V = [VMs(:,1:m)];
        else
            V = [VMs(:,1:m) MDs(:,1:n) DS(:,1:k)];
        end
        %V=[ MDs(:,1:6)];
        %mass normalization
        for ii = 1 : size(V, 2)
            V(:,ii) = V(:,ii) / (V(:,ii)'*(BeamAssembly.mass_matrix)*V(:,ii));
        end
        BeamReducedAssemblyDefected  = ReducedAssembly(myMesh,V);
        
        
        BeamReducedAssemblyDefected.DATA.M = BeamReducedAssemblyDefected.mass_matrix();
        BeamReducedAssemblyDefected.DATA.C = BeamReducedAssemblyDefected.damping_matrix(0,0,u0);
        BeamReducedAssemblyDefected.DATA.K =  BeamReducedAssemblyDefected.stiffness_matrix(u0);
        BeamReducedAssemblyDefected.DATA.Ud=U*xi;

        %% algorithm 3d
%         V_H=V(:,1:m)
%         Lin_sol=TI_lin.Solution.u;
%         Lin_sol_snap=[];
%         for ii=1:size(Lin_sol,2)
%             if rem(ii,10)==0
%                 Lin_sol_snap=[Lin_sol_snap Lin_sol(:,ii)];
%             end
%         end
%         eta=V_H'*Lin_sol_snap;
%         
%         Theta = QM_Theta_from_SMDs(V(:,m+1:m+n), MDs_names(1:n,:));
%         Xi = DS_Xi_QM(V(:,m+n+1:m+n+k), DSnames);
%         %QM uplifting
%         uu = zeros(size(V_H,1), size(eta, 2));
%         for tt = 1 : size(eta, 2)
%             uu(:,tt) = einsum('Iij,iJ,jK->IJK', Theta, eta(:,tt), eta(:,tt));
%         end
%         %
%                 uud = zeros(size(V_H,1), size(eta, 2));
%         for tt = 1 : size(eta, 2)
%             for ii= 1:size(xi,2)
%                 uud(:,tt) =uud(:,tt)+ einsum('Iij,iJ,jK->IJK', Xi, eta(:,tt), xi(:,ii));
%             end
%         end
%         
%         u_lin_ECSW = V_H * eta + 1/2*uu+uud;
%      
%         %Construct Gb
%         qq=(V'*V)^-1*V'*u_lin_ECSW;
%         tic
%         [G,b]=BeamReducedAssemblyDefected.constructGb(qq);
%         GBconstructTime=toc;
        %fNNLS
        tic
        [x_fnnls,w_fnnls]=fnnls(G'*G,G'*b,(norm(b)*0.01));
        fnnlsTime=toc
        tic;[x_nnls,w_nnls,info]=nnls(full(G),full(b),struct('Accy',1,'Tol',1/((norm(b)*kkk))));
        nnlTime=toc
        options = optimset('TolX',(norm(b)*0.01));
        %tic;x_lsq=lsqnonneg(G,b,options);
        lsqTime=toc
        nnz(x_fnnls)
        nnz(x_nnls)
        numofELEMENTS=[numofELEMENTS  nnz(x_nnls)];
        nnz(x_lsq)
        BeamReducedAssemblyDefected.DATA.elementWeights=x_nnls;
        %%
        f0_ROM_i= sort(sqrt(eig(BeamReducedAssemblyDefected.DATA.M\BeamReducedAssemblyDefected...
            .DATA.K))/2/pi) ;
        f0c_ROMHd=[f0c_ROMHd f0_ROM_i];
        %
        q0 = zeros(t,1);
        qd0 = zeros(t,1);
        qdd0 = zeros(t,1);
        
        
        TI_lin_red = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);
        
        % Modal linear Residual evaluation function handle
        Residual_lin_red = @(q,qd,qdd,t)residual_reduced_linear(q,qd,qdd,t,BeamReducedAssemblyDefected,F_ext);
        
        % time integration
        tic
        TI_lin_red.Integrate(q0,qd0,qdd0,tmax,Residual_lin_red);
        TI_lin_red.Solution.u = V * TI_lin_red.Solution.q;
        durationL=toc
        LdurationHd=[LdurationHd durationL];
        LtimeHd=[LtimeHd TI_lin_red.Solution.time'];
        LROMHd=[LROMHd TI_lin_red.Solution.u(dof,:)'];
        
        %% Reduced solution Noninear
        % For demonstration purposes, we simply reduce the nonlinear system using
        % out-of-plane bending modes. This is expected to produce bad results when
        % in-plane stretching is involved in the response.
        
        
        %TI_NL_alpha_red = GeneralizedAlpha('timestep',h,'rho_inf',0.7);
        TI_NL_alpha_red = ImplicitNewmark('timestep',h,'alpha',0.005);
        
        % Modal nonlinear Residual evaluation function handle
        Residual_NL_red = @(q,qd,qdd,t)residual_reduced_nonlinear_hyper_defected(q,qd,qdd,t,BeamReducedAssemblyDefected,F_ext);
        
        % time integration
        tic
        TI_NL_alpha_red.Integrate(q0,qd0,qdd0,tmax,Residual_NL_red);
        TI_NL_alpha_red.Solution.u = V * TI_NL_alpha_red.Solution.q;
        durationNL=toc
        
        NLdurationHd=[NLdurationHd durationNL];
        
        NLtimeHd=[NLtimeHd TI_NL_alpha_red.Solution.time'];
        NLROMHd=[NLROMHd TI_NL_alpha_red.Solution.u(dof,:)'];
        
        countHd=countHd+1;
        VbModeHd=[VbModeHd m];
        MDeriHd=[MDeriHd n];
          DSentiHd=[DSentiHd k];
          
          
          
          
    end
    end
end
[max_size, max_index] = max(cellfun('size', f0c_ROMHd, 1));
f0_ROM=zeros(max_size,size(f0c_ROMHd,2));
for ii=1:max_index
    f0c=f0c_ROMHd{ii};
    s=length(f0c);
    f0_ROM(1:s,ii)=f0c;
end

for i=1:countHd
    VmodeHd=VbModeHd(i);
    MderiHd=MDeriHd(i);
    DSsentiHd=DSentiHd(i)
      plot(LtimeHd(:,i), LROMHd(:,i),'DisplayName', ['hyperDPROM L (Newmark)  ' ' VM=' num2str(VmodeHd) ' MD=' num2str(MderiHd) ' DS=' num2str(DSsentiHd) ' time=' num2str(LdurationHd(i)) 's'])
    xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
    plot(NLtimeHd(:,i), NLROMHd(:,i),'DisplayName', ['HyperDPROM NL (Newmark)  ' ' VM=' num2str(VmodeHd) ' MD=' num2str(MderiHd) ' DS=' num2str(MderiHd) ' time=' num2str(NLdurationHd(i)) 's'])
    xlabel('time'); ylabel('u'); grid on; axis tight; legend('show')
end
end
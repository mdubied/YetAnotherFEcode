%to be deleted
%% POD plot
[UU SS VV]=svd(Lin_sol) ;
xSS=linspace(1,10002,10002);
figure
hold on
for ii=1:674
scatter(xSS(ii),SS(ii,ii))

end
for mod=1:2
    % mod = 4;
    elementPlot = elements(:,1:4); % plot only corners (otherwise it's a mess)
    figure('units','normalized','position',[.2 .1 .6 .8])
    
    PlotMesh(nodes, elementPlot, 0);
    v1 = reshape(UU(:,mod), 2, []).';
    PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor',Ly*1.1);
    title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0(mod),3) ' Hz'])
    
end
%%
%% Reduced solution
LROM=[];
Ltime=[];
Lduration=[];
NLROM=[];
NLtime=[];
NLduration=[];
f0c_ROM={};

nfsense = find_node(Lx/2,Ly/2,[],nodes); % node where to see the results
dof = get_index(nfsense, myMesh.nDOFPerNode );
dof=dof(2);

count=0;
VbMode=[];
MDeri=[];
for m=3;% m=[2 3 7 10]
    for n=6 % n=[0 m*(m+1)/2]
        
        t=m+n;
        % normalization
        % for ii = 1:m
        %     VMs(:,ii) = VMs(:,ii)/max(sqrt(sum(VMs(:,ii).^2,2)));
        % end
        % for ii = 1:n
        %    MDs(:,ii) = MDs(:,ii)/max(sqrt(sum(MDs(:,ii).^2,2)));
        % end
        % VMs = self.unconstrain_vector(VMs);
        
        if n==0
            V = [VMs(:,1:m)];
        else
            V = [VMs(:,1:m) MDs(:,1:n)];
        end
        %V=[ MDs(:,1:6)];
        %mass normalization
        for ii = 1 : size(V, 2)
            V(:,ii) = V(:,ii) / (V(:,ii)'*(BeamAssembly.mass_matrix)*V(:,ii));
        end
        BeamReducedAssembly  = ReducedAssembly(myMesh,V);
        
        
        BeamReducedAssembly.DATA.M = BeamReducedAssembly.mass_matrix();
        BeamReducedAssembly.DATA.C = BeamReducedAssembly.damping_matrix(0,0,u0);
        BeamReducedAssembly.DATA.K =  BeamReducedAssembly.stiffness_matrix(u0);
        BeamReducedAssembly.DATA.elementWeights=x_nnls;
        f0_ROM_i= sort(sqrt(eig(BeamReducedAssembly.DATA.M\BeamReducedAssembly...
            .DATA.K))/2/pi) ;
        f0c_ROM=[f0c_ROM f0_ROM_i];
        %
        q0 = zeros(t,1);
        qd0 = zeros(t,1);
        qdd0 = zeros(t,1);
        
        
        TI_lin_red = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);
        
        % Modal linear Residual evaluation function handle
        Residual_lin_red = @(q,qd,qdd,t)residual_reduced_linear(q,qd,qdd,t,BeamReducedAssembly,F_ext);
        
        % time integration
        tic
        TI_lin_red.Integrate(q0,qd0,qdd0,tmax,Residual_lin_red);
        TI_lin_red.Solution.u = V * TI_lin_red.Solution.q;
        durationL=toc
        Lduration=[Lduration durationL];
        Ltime=[Ltime TI_lin_red.Solution.time'];
        LROM=[LROM TI_lin_red.Solution.u(dof,:)'];
        
        %% Reduced solution Noninear
        % For demonstration purposes, we simply reduce the nonlinear system using
        % out-of-plane bending modes. This is expected to produce bad results when
        % in-plane stretching is involved in the response.
        
        
        %TI_NL_alpha_red = GeneralizedAlpha('timestep',h,'rho_inf',0.7);
        TI_NL_alpha_red = ImplicitNewmark('timestep',h,'alpha',0.005);
        
        % Modal nonlinear Residual evaluation function handle
        Residual_NL_red = @(q,qd,qdd,t)residual_reduced_nonlinear_hyper(q,qd,qdd,t,BeamReducedAssembly,F_ext);
        
        % time integration
        tic
        TI_NL_alpha_red.Integrate(q0,qd0,qdd0,tmax,Residual_NL_red);
        TI_NL_alpha_red.Solution.u = V * TI_NL_alpha_red.Solution.q;
        durationNL=toc
        
        NLduration=[NLduration durationNL];
        
        NLtime=[NLtime TI_NL_alpha_red.Solution.time'];
        NLROM=[NLROM TI_NL_alpha_red.Solution.u(dof,:)'];
        
        count=count+1;
        VbMode=[VbMode m];
        MDeri=[MDeri n];
    end
end
[max_size, max_index] = max(cellfun('size', f0c_ROM, 1));
f0_ROM=zeros(max_size,size(f0c_ROM,2));
for ii=1:max_index
    f0c=f0c_ROM{ii};
    s=length(f0c);
    f0_ROM(1:s,ii)=f0c;
end
%% Reduced solution defected
LROM=[];
Ltime=[];
Lduration=[];
NLROM=[];
NLtime=[];
NLduration=[];
f0c_ROM={};

nfsense = find_node(Lx/2,Ly/2,[],nodes); % node where to see the results
dof = get_index(nfsense, myMesh.nDOFPerNode );
dof=dof(2);

count=0;
VbMode=[];
MDeri=[];
for m=3;% m=[2 3 7 10]
    for n=6 % n=[0 m*(m+1)/2]
        for k=3
        
        t=m+n+k;
        % normalization
        % for ii = 1:m
        %     VMs(:,ii) = VMs(:,ii)/max(sqrt(sum(VMs(:,ii).^2,2)));
        % end
        % for ii = 1:n
        %    MDs(:,ii) = MDs(:,ii)/max(sqrt(sum(MDs(:,ii).^2,2)));
        % end
        % VMs = self.unconstrain_vector(VMs);
        
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
        BeamReducedAssembly  = ReducedAssembly(myMesh,V);
        
        
        BeamReducedAssembly.DATA.M = BeamReducedAssembly.mass_matrix();
        BeamReducedAssembly.DATA.C = BeamReducedAssembly.damping_matrix(0,0,u0);
        BeamReducedAssembly.DATA.K =  BeamReducedAssembly.stiffness_matrix(u0);
        BeamReducedAssembly.DATA.elementWeights=x_fnnls;
        f0_ROM_i= sort(sqrt(eig(BeamReducedAssembly.DATA.M\BeamReducedAssembly...
            .DATA.K))/2/pi) ;
        f0c_ROM=[f0c_ROM f0_ROM_i];
        %
        q0 = zeros(t,1);
        qd0 = zeros(t,1);
        qdd0 = zeros(t,1);
        
        
        TI_lin_red = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);
        
        % Modal linear Residual evaluation function handle
        Residual_lin_red = @(q,qd,qdd,t)residual_reduced_linear(q,qd,qdd,t,BeamReducedAssembly,F_ext);
        
        % time integration
        tic
        TI_lin_red.Integrate(q0,qd0,qdd0,tmax,Residual_lin_red);
        TI_lin_red.Solution.u = V * TI_lin_red.Solution.q;
        durationL=toc
        Lduration=[Lduration durationL];
        Ltime=[Ltime TI_lin_red.Solution.time'];
        LROM=[LROM TI_lin_red.Solution.u(dof,:)'];
        
        %% Reduced solution Noninear
        % For demonstration purposes, we simply reduce the nonlinear system using
        % out-of-plane bending modes. This is expected to produce bad results when
        % in-plane stretching is involved in the response.
        
        
        %TI_NL_alpha_red = GeneralizedAlpha('timestep',h,'rho_inf',0.7);
        TI_NL_alpha_red = ImplicitNewmark('timestep',h,'alpha',0.005);
        
        % Modal nonlinear Residual evaluation function handle
        Residual_NL_red = @(q,qd,qdd,t)residual_reduced_nonlinear_hyper(q,qd,qdd,t,BeamReducedAssembly,F_ext);
        
        % time integration
        tic
        TI_NL_alpha_red.Integrate(q0,qd0,qdd0,tmax,Residual_NL_red);
        TI_NL_alpha_red.Solution.u = V * TI_NL_alpha_red.Solution.q;
        durationNL=toc
        
        NLduration=[NLduration durationNL];
        
        NLtime=[NLtime TI_NL_alpha_red.Solution.time'];
        NLROM=[NLROM TI_NL_alpha_red.Solution.u(dof,:)'];
        
        count=count+1;
        VbMode=[VbMode m];
        MDeri=[MDeri n];
    end
    end
end
[max_size, max_index] = max(cellfun('size', f0c_ROM, 1));
f0_ROM=zeros(max_size,size(f0c_ROM,2));
for ii=1:max_index
    f0c=f0c_ROM{ii};
    s=length(f0c);
    f0_ROM(1:s,ii)=f0c;
end

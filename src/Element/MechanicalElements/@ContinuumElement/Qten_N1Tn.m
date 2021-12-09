% DpROM, reference: 
% Marconi et al. (2021). "A higher-order parametric nonlinear 
% reduced-order model for imperfect structures using Neumann 
% expansion". Nonlinear Dynamics. 
% https://doi.org/10.1007/s11071-021-06496-y

function Q = Qten_N1Tn(self, Ve, Ue, data)
    L1 = data.L1;
    L2 = data.L2;          
    C0 = self.initialization.C;
    H = self.initialization.H;
    Q = data.Qinit;

    X = self.quadrature.X;
    W = self.quadrature.W;
    for ii = 1:length(W)
        Xi = X(:, ii);  % quadrature points
        we = W(ii);     % quadrature weights
        [G, detJ] = shape_function_derivatives(self, Xi);               
        Y = G*Ue;
        G = G*Ve;
        C = C0*(we * detJ); % <--- includes we and detJ!

        % Q2n, nominal
        HG = H*G;
        Q2n_i = HG'*C*HG;

        % Q3dd, defected
        CHG = C*HG;
        L2YG = einsum('kla,aK,lJ->kJK',L2,Y,G);
        Q3d_i = einsum('kI,kJK->IJK',CHG,L2YG);
        Q3d_i = Q3d_i + permute(Q3d_i,[2 1 3]);

        % Q4dd, defected
        Q4dd_i = einsum('jIK,jk,kJL->IJKL',L2YG,C,L2YG);

        % Q3n, nominal 
        L1GG = einsum('kla,aK,lJ->kJK',L1,G,G);
        Q3n_i = einsum('kI,kJK->IJK',CHG,L1GG);
        Q3n_i = 0.5 * Q3n_i + permute(Q3n_i,[2 1 3]);

        % Q4d, defected
        Q4d_a = einsum('jIL,jk,kJK->IJKL',L2YG,C,L1GG);
        Q4d_i = 0.5*Q4d_a + permute(Q4d_a,[3 2 1 4]);

        % Q4n, nominal
        Q4n_i = 0.5 * einsum('jIK,jk,kLJ->IJKL',L1GG,C,L1GG);

        Q.Q2n{1}  = Q.Q2n{1}  + Q2n_i;
        Q.Q3d{1}  = Q.Q3d{1}  + Q3d_i;
        Q.Q4dd{1} = Q.Q4dd{1} + Q4dd_i;
        Q.Q3n{1}  = Q.Q3n{1}  + Q3n_i;
        Q.Q4d{1}  = Q.Q4d{1}  + Q4d_i;
        Q.Q4n{1}  = Q.Q4n{1}  + Q4n_i;
    end
end
        
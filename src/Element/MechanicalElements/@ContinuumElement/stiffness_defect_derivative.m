function dKdef = stiffness_defect_derivative(self, x, formulation)
% DpROM, reference: 
% Marconi et al. (2021). "A higher-order parametric nonlinear 
% reduced-order model for imperfect structures using Neumann 
% expansion". Nonlinear Dynamics. 
% https://doi.org/10.1007/s11071-021-06496-y
%
% tangent stiffness matrix derivative wrt the (scalar)
% amplitude xi_j of the defect vector Ud_j.
% dKdp is used to compute Defect Sensitivities.
% x --> defect
% formulation --> 'N0' or 'N1' (neumann expansion, order 0/1)

    u_defect = self.extract_element_data(x);
    X = self.quadrature.X;
    W = self.quadrature.W;            
    C = self.initialization.C;
    H = self.initialization.H;
    dKdef = self.initialization.K;
    switch upper(formulation)
        case 'N0'           % neumann order 0 (budiansky)
            A2fun = self.initialization.Afun;
        case {'N1','N1T'}	% neumann order 1/1t (t=truncated)
            A2fun = @(thd) A2_fun(self, thd);
    end
    for ii = 1:length(W)
        Xi = X(:, ii);  % quadrature points
        we = W(ii);     % quadrature weights
        [G,detJ] = shape_function_derivatives(self, Xi);
        thd = G*u_defect;
        A2 = A2fun(thd);
        int_dK = G'*(H'*C*A2 + A2'*C*H)*G; % Neumann formulation
        dKdef = dKdef + int_dK * detJ * we;
    end
end
        
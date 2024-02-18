function [ r, drdqdd,drdqd,drdq, c0] = residual_nonlinear_actu(q,qd,qdd,...
    t,Assembly,fActu,actuLeft,actuRight,actuSignalLeft,actuSignalRight)

% getting data 
M = Assembly.DATA.M;
C = Assembly.DATA.C;
u = Assembly.unconstrain_vector(q);
[K, F] = Assembly.tangent_stiffness_and_force(u);

% these matrices and the external forcing vector are appropriately constrained 
% according to the boundary conditions:
M_red = Assembly.constrain_matrix(M);
C_red = Assembly.constrain_matrix(C);
K_red = Assembly.constrain_matrix(K);

% forces
F_inertial = M_red * qdd;
F_damping = C_red * qd;
F_elastic = Assembly.constrain_vector(F);
F_external =  Assembly.constrain_vector(fActu(t,u));

% residual
r = F_inertial + F_damping + F_elastic - F_external ;

% residual derivatives
drdqdd = M_red;
drdqd = C_red;
drdq = K_red ...
    - actuSignalRight(t)*Assembly.constrain_matrix(actuRight.B2) ...
    - actuSignalLeft(t)*Assembly.constrain_matrix(actuLeft.B2);

% comparison norm of residual
c0 = norm(F_inertial) + norm(F_damping) + norm(F_elastic) + norm(F_external);
end
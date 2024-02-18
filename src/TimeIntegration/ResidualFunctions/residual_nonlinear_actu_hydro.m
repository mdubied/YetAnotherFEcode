function [ r, drdqdd,drdqd,drdq, c0] = residual_nonlinear_actu_hydro(q,qd,qdd,...
    t,Assembly,fActu,fTail,fSpine,fDrag,actuLeft,actuRight,actuSignalLeft,actuSignalRight)

% getting data 
M = Assembly.DATA.M;
C = Assembly.DATA.C;
u = Assembly.unconstrain_vector(q);
ud = Assembly.unconstrain_vector(qd);
udd = Assembly.unconstrain_vector(qdd);
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
F_external = Assembly.constrain_vector(fActu(t,u)) + ...
    Assembly.constrain_vector(fTail(u,ud)) + ...
    Assembly.constrain_vector(fSpine(u,ud,udd)) + ...
    Assembly.constrain_vector(fDrag(ud));

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
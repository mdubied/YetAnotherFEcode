function ydot = dynamics_for_RK4(y, t, Assembly, fActu, fTail, fSpine, fDrag)
    % Number of contrained degrees of freedom
    n_DOFs_constrained = length(y) / 2;   % Split y into position and velocity parts
    q = y(1:n_DOFs_constrained);          % Position vector, constrained
    qd = y(n_DOFs_constrained+1:end);     % Velocity vector, constrained

    % Get system matrices and unconstrain vectors for x, xd
    M = Assembly.DATA.M;
    C = Assembly.DATA.C;
    u = Assembly.unconstrain_vector(q);
    ud = Assembly.unconstrain_vector(qd);

    % Calculate stiffness and internal force vector
    [K, F] = Assembly.tangent_stiffness_and_force(u);
    
    % Constrain the matrices and vectors according to boundary conditions
    M_red = Assembly.constrain_matrix(M);
    C_red = Assembly.constrain_matrix(C);
    K_red = Assembly.constrain_matrix(K);
    F_elastic = Assembly.constrain_vector(F);

    % External force calculations
    F_actuation = Assembly.constrain_vector(fActu(t, u));
    F_tail = Assembly.constrain_vector(fTail(u, ud));
    F_spine = Assembly.constrain_vector(fSpine(u, ud, zeros(size(u)))); % TODO: change udd, or take it on the LHS of the equation
    F_drag = Assembly.constrain_vector(fDrag(ud));

    % Total external force
    F_external = F_actuation + F_tail + F_spine + F_drag;

    % Compute acceleration: xdd = M_red \ (F_external - F_damping - F_elastic)
    F_damping = C_red * qd;
    qdd = M_red \ (F_external - F_damping - F_elastic);

    % Convert second-order dynamics to first-order in state-space form
    ydot = [qd; qdd];
end

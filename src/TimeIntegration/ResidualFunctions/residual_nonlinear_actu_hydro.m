function [ r, drdqdd,drdqd,drdq, c0] = residual_nonlinear_actu_hydro(q,qd,qdd,...
    t,Assembly,fActu,fTail,fSpine,fDrag,actuLeft,actuRight,actuSignalLeft,actuSignalRight,fTailProp,fSpineProp,elements)
% ,fTailProp,fSpineProp,fDragProp,R,x0



% getting data 
M = Assembly.DATA.M;
C = Assembly.DATA.C;
u = Assembly.unconstrain_vector(q);
ud = Assembly.unconstrain_vector(qd);
udd = Assembly.unconstrain_vector(qdd);
[K, F] = Assembly.tangent_stiffness_and_force(u);

% for debugging (plot)
% nodes =  Assembly.Mesh.nodes;
% elementPlot = elements(:,1:4);
% u_plot = reshape(u, 3, []).';

% these matrices and the external forcing vector are appropriately constrained 
% according to the boundary conditions:
M_red = Assembly.constrain_matrix(M);
C_red = Assembly.constrain_matrix(C);
K_red = Assembly.constrain_matrix(K);

% forces
F_inertial = M_red * qdd;
F_damping = C_red * qd;
F_elastic = Assembly.constrain_vector(F);
F_external = Assembly.constrain_vector(fActu(t,u)) + ...  % factu as function of t only for the simple tip force
    Assembly.constrain_vector(fTail(u,ud)) + ...
    Assembly.constrain_vector(fSpine(u,ud,udd)) + ...
    Assembly.constrain_vector(fDrag(ud));

% fprintf('Max u: %.4f\n',max(u))
% fprintf('Min u: %.4f\n',min(u))

% fprintf('Max internal force: %.4f\n',norm(F_elastic))
% fprintf('Max actuation force: %.4f\n',norm(Assembly.constrain_vector(fActu(t,u))))
% % disp(max(actuSignalRight(t)*Assembly.constrain_matrix(actuRight.B2) - actuSignalLeft(t)*Assembly.constrain_matrix(actuLeft.B2)))

if norm(F_inertial + F_damping + F_elastic)>100
    fprintf('Warning')
end

% disp(norm(Assembly.constrain_vector(fSpine(u,ud,udd))))
% disp(norm(Assembly.constrain_vector(fActu(t,u))))

% residual
r = F_inertial + F_damping + F_elastic - F_external ;

% derivative actuation force
% der_actu = tip_actuation_force_derivatives_FOM(t,q,Assembly,fTailProp);

% derivatives tail force
der_tail_force = tail_force_derivatives_FOM(u,ud,Assembly,fTailProp);

% derivatives spine force
der_spine_force = spine_force_derivatives_FOM(u,ud,udd,fSpineProp.tensors);

% residual derivatives
drdqdd = M_red ;%-Assembly.constrain_matrix(der_spine_force.dfdqdd);
drdqd = C_red ;%- Assembly.constrain_matrix(der_spine_force.dfdqd)- Assembly.constrain_matrix(der_tail_force.dfdqd);% 
drdq = K_red ... %- actuSignalLeft*Assembly.constrain_matrix(der_actu.dfdq)...
     - actuSignalRight(t)*Assembly.constrain_matrix(actuRight.B2) ...
     - actuSignalLeft(t)*Assembly.constrain_matrix(actuLeft.B2) ;%...
     %- Assembly.constrain_matrix(der_spine_force.dfdq) ...
     %- Assembly.constrain_matrix(der_tail_force.dfdqd);% ...
     %- Assembly.constrain_matrix(der_spine_force.dfdq);
     
% fprintf('Max drdq: %.4f\n',normest(drdq));

% comparison norm of residual
c0 = norm(F_inertial) + norm(F_damping) + norm(F_elastic) + norm(F_external);
end
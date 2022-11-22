% residual_linear_sens
%
% Synthax:
% [r, drdsdd, drdsd, drds] = @(s,sd,sdd,t,it)residual_linear_sens(s,sd,sdd,t,ROMn_Assembly,qsol,qdsol,pd_fext,pd_fint,h)
%
% Description: as for the other residuals, this function is used as handle
% for the Newmark integration scheme. It is used to solve the ODE for the
% sensitivity, and therefore includes terms that need to be evaluated at
% each time step with the solution `qsol' and `qdsol'.
%
% INPUTS: 
% (1) s, sd, sdd, t:        variables for the function handle
% (2) ROMn_Assembly:        ROM-n assembly, needed to get the mass matrix M
% (3) qsol,qdsol,qddsol:    solutions of the simulations
% (4) pd_fext:              partial derivatives of external forces 
%                           (hydrodynamic forces) used for the sensitivity ODE
% (5) pd_fint:              partial derivative of internal forces
% (6) h:                    time step size h                       
%
% OUTPUTS:
% (1) [r, drdsdd, drdsd, drds]:     function handle describing the
%                                   residual and its partial derivatives
%     
% Additional notes:
%   The EoMs are
%   Mq_dd + Cq_d + K(p)q + f_int_nl(q,p) = f_ext(q,qd,p) % does K depends
%   on p?
%
% Last modified: 17/11/2022, Mathieu Dubied, ETH Zürich
function [r, drdsdd, drdsd, drds] = residual_linear_sens(s,sd,sdd,t,ROMn_Assembly,qsol,qdsol,qddsol,pd_fext,pd_fint,h)

    % compute current time step
    it = cast(t/h,"int16");

    % collect data
    M = ROMn_Assembly.DATA.M;
    C = ROMn_Assembly.DATA.C;
    K = ROMn_Assembly.DATA.K;
    qsolIt = qsol(:,it);
    qdsolIt = qdsol(:,it);
    qddsolIt = qddsol(:,it);
    derivative_fext_PROM = pd_fext(qsolIt,qdsolIt);
    derivative_fint_PROM = pd_fint(qsolIt);
    % external forces (hydrodynamic forces)
    dfextdq = derivative_fext_PROM.dfdq;
    dfextdqd = derivative_fext_PROM.dfdqd;
    dfextdp = derivative_fext_PROM.dfdp;
    % internal forces (only a function of q and not qd)
    dfintdq = derivative_fint_PROM.dfdq;
    dfintdp = derivative_fint_PROM.dfdp;
    dMdp = derivative_fint_PROM.dMdp;
    
    % residual, from sensitivity ODE. Partial derivative of fint are on the
    % RHS of the equation according to the formulation of previous papers
    r =  dMdp*qddsolIt + M*sdd - dfextdqd*sd - dfextdq*s  - dfextdp + ...
        dfintdq*s + dfintdp + C*sd + K*s;
    drdsdd = M;
    drdsd = -dfextdqd + C;
    drds = -dfextdq + dfintdq + K;
    
    
end
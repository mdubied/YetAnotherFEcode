% residual_linear_2ndSens
%
% Synthax:
% [r, drdsdd, drdsd, drds] = @(s,sd,sdd,t,it)residual_linear_2nSens(s,sd,sdd,t,ROMn_Assembly,qSol,sSol,pd_fext,pd_fint,h)
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
%   M(q)q_dd + Cq_d + Kq + f_int_nl(q,p) = f_ext(q,qd,p) 
%
% Last modified: 21/03/2023, Mathieu Dubied, ETH Zürich
function [r, drds2dd, drds2d, drds2] = residual_linear_2ndSens(s2,s2d,s2dd,t,ROMn_Assembly,qSol,sSol,pd_fext,pd_fint,h)

    % compute current time step
    it = cast(t/h,"int16");
    % number of shape variation parameters
    m = size(pd_fint(qSol.q(:,1)).dfdp,2);

    % collect data
    M = tensor(ROMn_Assembly.DATA.M);
    C = tensor(ROMn_Assembly.DATA.C);
    K = tensor(ROMn_Assembly.DATA.K);
    qsol = qSol.q;
    qdsol = qSol.qd;
    qddsol = qSol.qdd;
    ssol = sSol.q;
    sdsol = sSol.qd;
    sddsol = sSol.qdd;
    qsolIt = qsol(:,it);
    qdsolIt = qdsol(:,it);
    qddsolIt = qddsol(:,it);
    ssolIt = tensor(ssol(:,it));
    sdsolIt = tensor(sdsol(:,it));
    sddsolIt = tensor(sddsol(:,it));
    derivative_fext_PROM = pd_fext(qsolIt,qdsolIt);
    derivative_fint_PROM = pd_fint(qsolIt);

    % external forces (hydrodynamic forces)
    % first order derivatives

    dfextdq = tensor(derivative_fext_PROM.dfdq);

    
    dfextdqd = tensor(derivative_fext_PROM.dfdqd);
    dfextdp = derivative_fext_PROM.dfdp;
    % second order derivatives
    dfextdqdp = tensor(derivative_fext_PROM.dfdqdp);
    dfextdqddp = tensor(derivative_fext_PROM.dfdqddp);
    dfextdqdq = tensor(derivative_fext_PROM.dfdqdq);
    dfextdqdqd = tensor(derivative_fext_PROM.dfdqdqd);
    dfextdqddq = tensor(derivative_fext_PROM.dfdqddq);
    dfextdqddqd = tensor(derivative_fext_PROM.dfdqddqd);
    dfextdpdq = tensor(derivative_fext_PROM.dfdpdq);
    dfextdpdqd = tensor(derivative_fext_PROM.dfdpdqd);


    % internal forces (only a function of q and not qd)
    
    dfintdp = tensor(derivative_fint_PROM.dfdp);
    
    if m == 1
        dfintdq = tensor(derivative_fint_PROM.dfdq,size(derivative_fint_PROM.dfdq));
        dfintdqdp = tensor(derivative_fint_PROM.dfdqdp,[size(derivative_fint_PROM.dfdqdp) 1]);
        dMdp = tensor(derivative_fint_PROM.dMdp,[size(derivative_fint_PROM.dMdp) 1]);
        dfintdpdp = tensor(derivative_fint_PROM.dfdp2,[size(derivative_fint_PROM.dfdp2) 1]);
        dfintdpdq = tensor(derivative_fint_PROM.dfdpdq,[size(derivative_fint_PROM.dfdpdq) 1]);
    else
        dfintdq = tensor(derivative_fint_PROM.dfdq);
        dfintdqdp = tensor(derivative_fint_PROM.dfdqdp);
        dMdp = tensor(derivative_fint_PROM.dMdp);
        dfintdpdp = tensor(derivative_fint_PROM.dfdp2);
        dfintdpdq = tensor(derivative_fint_PROM.dfdpdq);
    end

    dfintdqdq = tensor(derivative_fint_PROM.dfdq2);
    
    
   
    

    % convert  s2d, s2dd to tensors
    if m == 1
        s2 = tensor(s2,[size(s2) 1]);
        s2d = tensor(s2d,[size(s2d) 1]);
        s2dd = tensor(s2dd,[size(s2dd) 1]);
    else
        s2 = tensor(s2);
        s2d = tensor(s2d);
        s2dd = tensor(s2dd);
    end
   
    
    % residual, from sensitivity ODE. Partial derivative of fint are on the
    % RHS of the equation according to the formulation of previous papers
    r =  double(2*ttt(dMdp,sddsolIt,2,1) + ttt(M,s2dd,2,1) + ttt(C,s2d,2,1)+ ...
            ttt(dfintdqdp,ssolIt,2,1) + ttt(ttt(dfintdqdq,ssolIt,3,1),ssolIt,2,1)  + ...
            ttt(dfintdq,s2,2,1) + dfintdpdp + ttt(dfintdpdq,ssolIt,2,1) - ...
            ttt(dfextdqdp,ssolIt,2,1) - ttt(ttt(dfextdqdq,ssolIt,3,1),ssolIt,2,1) - ...
            ttt(ttt(dfextdqdqd,sdsolIt,3,1),ssolIt,2,1) - ttt(dfextdq,s2,2,1) - ...
            ttt(dfextdqddp,sdsolIt,2,1) - ttt(ttt(dfextdqddq,ssolIt,3,1),sdsolIt,2,1) - ...
            ttt(ttt(dfextdqddqd,sdsolIt,3,1),sdsolIt,2,1) - ttt(dfextdqd,s2d,2,1) - ...
            ttt(dfextdpdq,ssolIt,3,1)- ttt(dfextdpdqd,sdsolIt,3,1) + ttt(K,s2,2,1));
    drds2dd = double(M);
    drds2d = double(-dfextdqd) + double(C);
    drds2 = double(-dfextdq) + double(dfintdq) + double(K);
    
    
end
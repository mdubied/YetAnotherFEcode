% residual_linear_sens
%
% Synthax:
% [r, drdsdd, drdsd, drds] = @(s,sd,sdd,t,it)residual_linear_sens(s,sd,sdd,t,ROMn_Assembly,qsol,qdsol,pd_fext,pd_fint,h,varargin)
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
% (4) pd_fhydro:            partial derivatives of used for the sensitivity
%                           ODE
% (5) pd_fint:              partial derivative of internal forces
% (6) h:                    time step size h 
% (7) varargin:             additional name-value pairs if actuation is
%                           present. See below for the parser function.
%
% OUTPUTS:
% (1) [r, drdsdd, drdsd, drds]:     function handle describing the
%                                   residual and its partial derivatives
%     
%
% Last modified: 29/10/2023, Mathieu Dubied, ETH Zürich
function [r, drdsdd, drdsd, drds] = ...
residual_linear_sens(s,sd,sdd,t,ROMn_Assembly,qsol,qdsol,qddsol, ...
        pd_fint,pd_tail,pd_spine,pd_drag,pd_actuTop,pd_actuBottom, ...
        actuSignalTop,actuSignalBottom,h)
    
    % compute current time step
    it = cast(t/h,"int16");

    % number of shape variation parameters
    m = size(pd_fint(qsol(:,1)).dfdp,2);

    % COLLECT DATA ________________________________________________________
    M = ROMn_Assembly.DATA.M;
    C = ROMn_Assembly.DATA.C;
    K = ROMn_Assembly.DATA.K;
    
    qsolIt = qsol(:,it);
    qdsolIt = qdsol(:,it);
    qddsolIt = qddsol(:,it);

    % EVALUATE FUNCTION HANDLE ____________________________________________
    aTop = actuSignalTop(t);
    aBottom = actuSignalBottom(t);
    
    der_fint = pd_fint(qsolIt);
    der_tail = pd_tail(qsolIt,qdsolIt);
    der_spine = pd_spine(qsolIt,qdsolIt,qddsolIt);
    der_drag = pd_drag(qdsolIt);
    der_actuTop = pd_actuTop(aTop);
    der_actuBottom = pd_actuBottom(aBottom);

    % GATHER PARTIAL DERIVATIVES __________________________________________
    % internal forces (only a function of q and not qd)
    dfintdq = der_fint.dfdq;
    dfintdp = der_fint.dfdp;
    if m == 1
        m1 = size(der_fint.dMdp,1);
        m2 = size(der_fint.dMdp,2);
        dMdp = tenzeros(m1,m2,1);
        dMdp(:,:,1) = der_fint.dMdp;
    else
        dMdp = tensor(der_fint.dMdp);
    end

    % tail pressure force
    dfTaildq = der_tail.dfdq;
    dfTaildqd = der_tail.dfdqd;
    dfTaildp = der_tail.dfdp;

    % spine change in momentum
    dfSpinedq = der_spine.dfdq;
    dfSpinedqd = der_spine.dfdqd;
    dfSpinedqdd = der_spine.dfdqdd;
    dfSpinedp = der_spine.dfdp;

    % drag force
    dfDragdqd = der_drag.dfdqd;
    dfDragdp = der_drag.dfdp;

    % actuation forces
    dfactTopdq = der_actuTop.dfdq;
    dfactTopdp = der_actuTop.dfdp;
    dfactBottomdq = der_actuBottom.dfdq;
    dfactBottomdp = der_actuBottom.dfdp;

    % COMPUTE RESIDUAL ____________________________________________________
    r =  double(ttv(tensor(dMdp),qddsolIt,2)) + M*sdd + C*sd + K*s + dfintdq*s + dfintdp ...
        - dfTaildq*s - dfTaildqd*sd -dfTaildp ...
        - dfSpinedq*s -dfSpinedqd*sd - dfSpinedqdd*sdd -dfSpinedp ...
        - dfDragdqd*sd - dfDragdp ...
        - dfactTopdq*s -dfactTopdp - dfactBottomdq*s - dfactBottomdp;
   
    drdsdd = M - dfSpinedqdd;
    drdsd = C - dfTaildqd - dfSpinedqd - dfDragdqd;
    drds = K + dfintdq - dfTaildq - dfSpinedq - dfactTopdq - dfactBottomdq ;
      
end

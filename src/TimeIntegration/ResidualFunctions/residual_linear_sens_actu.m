% residual_linear_sens_actu
%
% Synthax:
% [r, drdsdd, drdsd, drds] = @(s,sd,sdd,t,it)residual_linear_sens_actu(s,sd,sdd,t,ROMn_Assembly,qsol,qdsol,pd_fext,pd_fint,h,varargin)
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
% Last modified: 18/05/2024, Mathieu Dubied, ETH Zürich
function [r, drdsdd, drdsd, drds] = ...
residual_linear_sens_actu(s,sd,sdd,t,ROMn_Assembly,qsol,qdsol,qddsol, ...
        pd_actu,h)
    
    % compute current time step
    it = cast(t/h,"int16");

    % COLLECT DATA ________________________________________________________
    M = ROMn_Assembly.DATA.M;
    C = ROMn_Assembly.DATA.C;
    K = ROMn_Assembly.DATA.K;
    
    qsolIt = qsol(:,it);
    qdsolIt = qdsol(:,it);
    qddsolIt = qddsol(:,it);

    % EVALUATE FUNCTION HANDLE ____________________________________________
%     aTop = actuSignalTop(t);
%     aBottom = actuSignalBottom(t);
%     
%     der_actuTop = pd_actuTop(aTop);
%     der_actuBottom = pd_actuBottom(aBottom);
    der_actu = pd_actu(t,qsolIt);

    % GATHER PARTIAL DERIVATIVES __________________________________________
    % actuation forces
%     dfactTopdq = der_actuTop.dfdq;
%     dfactTopdp = der_actuTop.dfdp;
%     dfactBottomdq = der_actuBottom.dfdq;
%     dfactBottomdp = der_actuBottom.dfdp;
    dfactdq = der_actu.dfdq;
    dfactdp = der_actu.dfdp;

    % COMPUTE RESIDUAL ____________________________________________________
    r =  M*sdd + C*sd + K*s ...
        - dfactdq*s - dfactdp;
   
    drdsdd = M;
    drdsd = C;
    drds = K - dfactdq;
    
end

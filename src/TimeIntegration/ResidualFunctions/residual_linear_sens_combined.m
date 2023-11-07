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
% (1) r:                    function handle describing the residual
%
% Last modified: 07/11/2023, Mathieu Dubied, ETH Zürich
function r = residual_linear_sens_combined(s,sd,sdd,t, ...
        qsol,qdsol,qddsol, drdqdd, drdqd, drdq, ...
        pd_fint,pd_tail,pd_spine,pd_drag,pd_actuTop,pd_actuBottom, ...
        actuSignalTop,actuSignalBottom)
    

    % number of shape variation parameters
    m = size(pd_fint(qsol(:,1)).dfdp,2);

    % COLLECT DATA ________________________________________________________   
    qsolIt = qsol;
    qdsolIt = qdsol;
    qddsolIt = qddsol;

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

    dfTaildp = der_tail.dfdp;

    % spine change in momentum
    dfSpinedp = der_spine.dfdp;

    % drag force
    dfDragdp = der_drag.dfdp;

    % actuation forces
    dfactTopdp = der_actuTop.dfdp;
    dfactBottomdp = der_actuBottom.dfdp;


    % COMPUTE RESIDUAL ____________________________________________________
    r =  drdqdd*sdd + drdqd*sd + drdq*s ...
        + double(ttv(tensor(dMdp),qddsolIt,2)) +  dfintdp ...
        -dfTaildp - dfSpinedp - dfDragdp ...
        -dfactTopdp - dfactBottomdp;
       
end

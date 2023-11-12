% residual_linear_sens
%
% Synthax:
% r = residual_linear_sens_combined(s,sd,sdd,t, ...
%       qsol,qdsol,qddsol, drdqdd, drdqd, drdq, ...
%       pd_fint,pd_tail,pd_spine,pd_drag,pd_actuTop,pd_actuBottom, ...
%       actuSignalTop,actuSignalBottom)
%
% Description: this residual, passed as a handle to
% solve_EoMs_and_sensitivties, is used to compute the sensitivities in a
% combined fashion at each integration step of the EoMs. Details of the
% algorithm can be found in "Sensitivity analysis for dynamic mechanical systems
% with finite rotations" by O. Brüls and P. Eberhard.
%
% INPUTS: 
% (1) s, sd, sdd, t:        variables for the function handle
% (2) qsol,qdsol,qddsol:    solutions of the EoMs solver
% (3) drdqdd,drqqd,drdq:    residual of the EoMs at the considered time
%                           step
% (4) pd_fint:              partial derivatives of internal forces
% (5) pd_tail:              partial derivative of tail pressure force
% (6) pd_spine:             partial derivative of spine change in momentum
% (7) od_drag:              partial derivative of drag forces
% (8) pd_actuTop:           partial derivative of actuation (top/left muscle)
% (9) pd_actuBottom:        partial derivative of actuation (bottom/right muscle)
% (10) actuSignalTop:       actuation signal of top/left muscle as a
%                           function handle
% (12) actuSignalBottom:    actuation signal of bottom/right muscle as a
%                           function handle
%
% OUTPUTS:
% (1) r:                    function handle describing the residual
%
% Last modified: 12/11/2023, Mathieu Dubied, ETH Zürich
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

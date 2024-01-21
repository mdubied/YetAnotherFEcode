% residual_linear_sens_co
%
% Synthax:
% r = residual_linear_sens_combined_co(s,sd,sdd,t, ...
%       qsol,qdsol,qddsol, drdqdd, drdqd, drdq, ...
%       pd_fint,pd_tail,pd_spine,pd_drag,pd_actuTop,pd_actuBottom, ...
%       actuSignalTop,actuSignalBottom)
%
% Description: this residual, passed as a handle to
% solve_EoMs_and_sensitivties, is used to compute the sensitivities in a
% combined fashion at each integration step of the EoMs. Details of the
% algorithm can be found in "Sensitivity analysis for dynamic mechanical 
% systems with finite rotations" by O. Brüls and P. Eberhard.
% This version is used for co-optimisation problems (shape + actuation).
%
% INPUTS: 
% (1) s, sd, sdd, t:        variables for the function handle
% (2) qsol,qdsol,qddsol:    solutions of the EoMs solver
% (3) drdqdd,drqqd,drdq:    residual of the EoMs at the considered time
%                           step
% (4) pd_fint:              partial derivatives of internal forces
% (5) pd_tail:              partial derivative of tail pressure force
% (6) pd_spine:             partial derivative of spine change in momentum
% (7) pd_drag:              partial derivative of drag forces
% (8) pd_actu:              partial derivatives of the actuation force
% (9) pActu:                current values of the actuation parameters

%
% OUTPUTS:
% (1) r:                    function handle describing the residual
%
% Last modified: 07/01/2024, Mathieu Dubied, ETH Zürich
function r = residual_linear_sens_combined_co(s,sd,sdd,t, ...
        qsol,qdsol,qddsol, drdqdd, drdqd, drdq, ...
        pd_fint,pd_tail,pd_spine,pd_drag, ...
        pd_shape_actuTop,pd_shape_actuBottom, ...
        actuSignalTop,actuSignalBottom, ...
        pd_actu_actu,pActu)
    

    % number of shape variation parameters
    mShape = size(pd_fint(qsol(:,1)).dfdp,2);
    
    % number of actuation signal parameters
    mActu = length(pActu); 
    
    % total number of parameters (shape + actuation)
    mdTot = mShape + mActu;
    
    % size of reduced order basis
    m = size(pd_fint(qsol(:,1)).dfdp,1);

    % COLLECT DATA ________________________________________________________   
    qsolIt = qsol;
    qdsolIt = qdsol;
    qddsolIt = qddsol;

    % EVALUATE FUNCTION HANDLE ____________________________________________
    aTop = actuSignalTop(t);
    aBottom = actuSignalBottom(t);
    
    % partial derivatives w.r.t. shape
    der_fint = pd_fint(qsolIt);
    der_tail = pd_tail(qsolIt,qdsolIt);
    der_spine = pd_spine(qsolIt,qdsolIt,qddsolIt);
    der_drag = pd_drag(qdsolIt);
    der_shape_actuTop = pd_shape_actuTop(aTop);
    der_shape_actuBottom = pd_shape_actuBottom(aBottom);
    
    % augment partial derivatives with actuation dimensions
    der_fint.dfdp = [der_fint.dfdp,zeros(m,mActu)];
    der_fint.dMdp(:,:,mShape+1:mdTot) = zeros(m,m,mActu);
    der_tail.dfdp = [der_tail.dfdp,zeros(m,mActu)];
    der_spine.dfdp = [der_spine.dfdp,zeros(m,mActu)];
    der_drag.dfdp = [der_drag.dfdp,zeros(m,mActu)];
    der_shape_actuTop.dfdp = [der_shape_actuTop.dfdp,zeros(m,mActu)];
    der_shape_actuBottom.dfdp = [der_shape_actuBottom.dfdp,zeros(m,mActu)];
    
    % partial derivative w.r.t. actuation
    der_actu_actu = pd_actu_actu(t,qsol,pActu);
    
    % augment partial derivative with shape dimension
    der_actu_actu.dfdp = [zeros(m,mShape),der_actu_actu.dfdp];

    % GATHER PARTIAL DERIVATIVES __________________________________________
    % internal forces (only a function of q and not qd)
    dfintdp = der_fint.dfdp;
    if mShape == 1
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
    dfactTopdp = der_shape_actuTop.dfdp;
    dfactBottomdp = der_shape_actuBottom.dfdp;
    dfactdpActu = der_actu_actu.dfdp;

    % COMPUTE RESIDUAL ____________________________________________________
    r =  drdqdd*sdd  + drdqd*sd + drdq*s ...
        + double(ttv(tensor(dMdp),qddsolIt,2)) +  dfintdp ...
        - dfTaildp - dfSpinedp - dfDragdp ...
        - dfactTopdp - dfactBottomdp - dfactdpActu;
       
end

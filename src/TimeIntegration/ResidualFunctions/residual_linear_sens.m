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
% Last modified: 28/03/2023, Mathieu Dubied, ETH Zürich
function [r, drdsdd, drdsd, drds] = residual_linear_sens(s,sd,sdd,t,ROMn_Assembly,qsol,qdsol,qddsol,pd_fhydro,pd_fint,h,varargin)
    
    % parse input
    [ACTUATION,pd_factTop,pd_factBottom] = parse_inputs(varargin{:});

    % compute current time step
    it = cast(t/h,"int16");

    % number of shape variation parameters
    m = size(pd_fint(qsol(:,1)).dfdp,2);

    % collect data
    M = ROMn_Assembly.DATA.M;
    C = ROMn_Assembly.DATA.C;
    K = ROMn_Assembly.DATA.K;
    qsolIt = qsol(:,it);
    qdsolIt = qdsol(:,it);
    qddsolIt = qddsol(:,it);
    derivative_fhydro_PROM = pd_fhydro(qsolIt,qdsolIt);
    derivative_fint_PROM = pd_fint(qsolIt);
    % hydrodynamic forces
    dfhydrodq = derivative_fhydro_PROM.dfdq;
    dfhydrodqd = derivative_fhydro_PROM.dfdqd;
    dfhdrodp = derivative_fhydro_PROM.dfdp;
    % internal forces (only a function of q and not qd)
    dfintdq = derivative_fint_PROM.dfdq;
    dfintdp = derivative_fint_PROM.dfdp;
    if m == 1
        dMdp = derivative_fint_PROM.dMdp;
    else
        dMdp = tensor(derivative_fint_PROM.dMdp);
    end
    % actuation forces
    if ACTUATION
        a1=1+0.004*sin(t*2*pi/5);
        a2=1-0.004*sin(t*2*pi/5);
        derivative_factTop_PROM = pd_factTop(a1);
        derivative_factBottom_PROM = pd_factBottom(a2);
        
        dfactTopdq = derivative_factTop_PROM.dfdq;
        dfactTopdp = derivative_factTop_PROM.dfdp;
        dfactBottomdq = derivative_factBottom_PROM.dfdq;
        dfactBottomdp = derivative_factBottom_PROM.dfdp;
    end


    % residual, from sensitivity ODE. 
    % without actuation
    if m == 1
        r =  dMdp*qddsolIt + M*sdd - dfhydrodqd*sd - dfhydrodq*s  - dfhdrodp + ...
            dfintdq*s + dfintdp + C*sd + K*s;
    else
        r =  double(ttv(dMdp,qddsolIt,2)) + M*sdd - dfhydrodqd*sd - dfhydrodq*s  - dfhdrodp + ...
            dfintdq*s + dfintdp + C*sd + K*s;
    end
    drdsdd = M;
    drdsd = -dfhydrodqd + C;
    drds = -dfhydrodq + dfintdq + K;
    % with actuation
    if ACTUATION
        r = r - dfactTopdq*s - dfactBottomdq*s- dfactTopdp - dfactBottomdp;
        drds = drds - dfactTopdq - dfactBottomdq;
    end
    
end

% parse input
function [ACTUATION,pd_factTop,pd_factBottom] = parse_inputs(varargin)
defaultACTUATION = 0;
defaultPd_factTop = 0;
defaultPd_factBottom = 0;
p = inputParser;

addParameter(p,'ACTUATION',defaultACTUATION,@(x)validateattributes(x, ...
                {'numeric'},{'nonempty','integer','positive'}) );
addParameter(p,'pd_factTop',defaultPd_factTop);
addParameter(p,'pd_factBottom',defaultPd_factBottom);

parse(p,varargin{:});

ACTUATION = p.Results.ACTUATION;
pd_factTop = p.Results.pd_factTop;
pd_factBottom = p.Results.pd_factBottom;

end
% residual_linear_sens_combined_actu
%
% Synthax:
% r = residual_linear_sens_combined_actu(s,sd,sdd,t,qsol,drdqdd,drdqd,drdq,pd_actu)
%
% Description: as for the other residuals, this function is used as handle
% for the Newmark integration scheme. It is used to solve the ODE for the
% sensitivity, and therefore includes terms that need to be evaluated at
% each time step with the solution `qsol'.
%
% INPUTS: 
% (1) s, sd, sdd, t:        variables for the function handle
% (2) qsol:                 solutions of the simulations
% (3) drdqdd,drdqd,drdq:    derivatives of the EoMs residual
% (4) pd_actu:              partial derivatives of the actuation force
%
% OUTPUTS:
% (1) r:                    function handle describing the residual
%
% Last modified: 10/11/2023, Mathieu Dubied, ETH Zürich
function r = residual_linear_sens_combined_actu(s,sd,sdd,t, ...
        qsol,drdqdd,drdqd,drdq,pd_actu)
    
    % EVALUATE FUNCTION HANDLE ____________________________________________   
    der_actu = pd_actu(t,qsol);
       
    % GATHER PARTIAL DERIVATIVES __________________________________________
    dfactdp = der_actu.dfdp;

    % COMPUTE RESIDUAL ____________________________________________________
    r =  drdqdd*sdd + drdqd*sd + drdq*s - dfactdp;
       
end

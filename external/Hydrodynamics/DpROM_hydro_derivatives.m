% DpROM_hydro_derivatives
%
% Synthax:
% der = DpROM_hydro_derivatives(q,qd,xi,tensors_DpROM)
%
% Description: This function returns the partial derivatives of the
% sensitity ODE, evaluated at the solution displacements q and velocities qdot.
%
% INPUTS: 
% (1) q:                vector of time domain displacements
% (2) qdot:             vector of time domain velocities
% (3) tensors_DpROM:    structure array containing reduced tensors of
%                       F_hydro model
%
% OUTPUTS:
% (1) der:              strucure array containing first order partial
%                       derivatives needed to solve the sensitivity ODE. 
%                       These derivatives are stored in the following 
%                       fields:'dfdp','dfdq',dfdqd'
%     
%
% Additional notes:
%   - q and qd should be understood as eta and dot{eta}.
%   - p should be understood as xi.
%
% Last modified: 17/11/2022, Mathieu Dubied, ETH Zürich
function der = DpROM_hydro_derivatives(q,qd,xi,tensors_DpROM,FOURTHORDER)

% populate tensors
Tr1 = tensors_DpROM.Tr1;
Tr2 = tensors_DpROM.Tr2; 
Tru2 = tensors_DpROM.Tru2;
Tru3 = tensors_DpROM.Tru3;
Trudot2 = tensors_DpROM.Trudot2;
Trudot3 = tensors_DpROM.Trudot3;
Truu3 = tensors_DpROM.Truu3;
Truudot3 = tensors_DpROM.Truudot3;
Trudotudot3 = tensors_DpROM.Trudotudot3;
if FOURTHORDER
    Truu4 = tensors_DpROM.Truu4;
    Truudot4 = tensors_DpROM.Truudot4;
    Trudotudot4 = tensors_DpROM.Trudotudot4;
end

% evaluate partial derivatives
df1dp = Tr2;
df2dp = ttv(Tru3,q,2) + ttv(Trudot3,qd,2);

dfdp = df1dp + df2dp;
if FOURTHORDER
    df3dp = einsum('IijL,i,j->IL',Truu4,q,q) + ...
            einsum('IijL,i,j->IL',Truudot4,q,qd) + ...
            einsum('IijL,i,j->IL',Trudotudot4,qd,qd);
    dfdp = dfdp + df3dp;
end

df1dq = zeros(size(Tr2,1));
df1dqd = zeros(size(Tr2,1));
df2dq = Tru2 + ttv(Tru3,xi,3);
df2dqd = Trudot2 + ttv(Trudot3,xi,3);
df3dq = ttv(Truu3,q,3) + ttv(Truu3,q,2) + ttv(Truudot3,qd,3);
df3dqd = ttv(Trudotudot3,qd,3) + ttv(Trudotudot3,qd,2) + ttv(Truudot3,q,2);
if FOURTHORDER
    df3dq = df3dq + einsum('ILij,i,j->IL',Truu4,q,xi) + ...
                    einsum('IiLj,i,j->IL',Truu4,q,xi) + ...
                    einsum('ILij,i,j->IL',Truudot4,qd,xi);
    df3dqd = df3dqd + einsum('ILij,i,j->IL',Trudotudot4,qd,xi) + ...
                    einsum('IiLj,i,j->IL',Trudotudot4,qd,xi) + ...
                    einsum('IiLj,i,j->IL',Truudot4,q,xi);
end


dfdq = df1dq + df2dq + df3dq;
dfdqd = df1dqd + df2dqd + df3dqd;

% store results in output struct
der.dfdp = double(dfdp);
der.dfdq = double(dfdq);
der.dfdqd = double(dfdqd);


end

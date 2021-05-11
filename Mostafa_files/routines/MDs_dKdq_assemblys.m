
function dKdq = MDs_dKdq_assemblys(model,Phi_i,n)

Nel  = model.Nelem;
str = pad(sprintf(' Assembling %d elements ...',Nel),33,'.');
fprintf(str)
t0=tic;

prop = model.properties.prop;
E = prop(1);                % elastic modulus
v = prop(2);                % poisson coefficient
coeff = E/((1+v)*(1-2*v));  % stiffness matrix constant
C = [1-v v v 0 0 0;         % constitutive law matrix
    v 1-v v 0 0 0;
    v v 1-v 0 0 0;
    0 0 0 .5*(1-2*v) 0 0;
    0 0 0 0 .5*(1-2*v) 0;
    0 0 0 0 0 .5*(1-2*v)]*coeff;

H = [1 0 0 0 0 0 0 0 0;
     0 0 0 0 1 0 0 0 0;
     0 0 0 0 0 0 0 0 1;
     0 1 0 1 0 0 0 0 0;
     0 0 1 0 0 0 1 0 0;
     0 0 0 0 0 1 0 1 0];

conn     = model.connectivity;
elements = model.elements;
nodes    = model.nodes;

nnel = size(elements,2)-1;  % number of nodes per element
neldofs = size(conn,2)-1;	% number of dofs per element
 
% Gauss quadrature points and weights
switch nnel
    case 4
        [X,W] = inttet(1);
    case 10
        [X,W] = inttet(n);
    case 20
        [x,w]=lgwt(n,-1,1);
        X = zeros(3,n^3);
        W = zeros(n^3,1);
        cont = 1;
        for ii = 1:n
            for jj = 1:n
                for kk = 1:n
                    X(:,cont) = [x(ii) x(jj) x(kk)].';
                    W(cont) = w(ii)*w(jj)*w(kk);
                    cont = cont+1;
                end
            end
        end
    otherwise
        error([ num2str(nnel) ' nodes per element found'])
end

% SMART ASSEMBLY ******************************************************** %
% Author: Daniele Giannini
% Last modified: 1-Oct-2018, by Jacopo Marconi

edofMat = conn(:,2:end);
iK = reshape(kron(edofMat,ones(neldofs,1))',neldofs^2*Nel,1);
jK = reshape(kron(edofMat,ones(1,neldofs))',neldofs^2*Nel,1);

ssK=zeros(neldofs^2,Nel);
for ee=1:Nel
    
    el_nodes = elements(ee,2:end);
    el_dofs  = conn(ee,2:end);
    
    xyz = nodes(el_nodes,2:4);
    Phi_i_el = Phi_i(el_dofs);
    
    dKdq_e = MDs_dKdq_elem(xyz,Phi_i_el,X,W,H,C);
    ssK(:,ee)=dKdq_e(:);
    
end
sK = reshape(ssK,neldofs^2*Nel,1);

dKdq = sparse(iK,jK,sK); 

fprintf(' %.2f s\n',toc(t0))
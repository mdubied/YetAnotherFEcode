
function dKdxi = stiffness_matrix_sensitivity(myAssembly, elements, U, formulation)

t0=tic;

% data from myAssembly
nodes    = myAssembly.Mesh.nodes;                   % nodes table
nel      = myAssembly.Mesh.nElements;           	% number of elements
nnodes   = myAssembly.Mesh.nNodes;               	% number of nodes
freedofs = myAssembly.Mesh.EBC.unconstrainedDOFs;   % free DOFs

% Element data (assuming ALL the elements have the same properties in terms
% of material and quadrature rule)
Element = myAssembly.Mesh.Elements(1).Object;   % first element
XGauss = Element.quadrature.X;                  % quadrature points
WGauss = Element.quadrature.W;                  % quadrature weights
C = Element.initialization.C;                   % constitutive matrix

DIM = size(nodes,2);                            % 2D/3D problem
ndofs = nnodes * DIM;                           % total number of DOFs
nfree = length(freedofs);                       % number of free DOFs

% build connectivity table
conn = zeros(nel, length(Element.iDOFs));
for ee = 1 : nel
    conn(ee, :) = myAssembly.Mesh.Elements(ee).Object.iDOFs;
end

elements = int64(elements); % convert into integers
conn = int64(conn);         % convert into integers

% INPUT SIZE CHECK ________________________________________________________
% for tensor assembly, always use the full DOFs displacement vectors
if size(U,1)==nfree
    u = zeros(ndofs,size(U,2));
    u(freedofs,:) = U;
    U = u;
elseif size(U,1)~=ndofs
    error('Wrong dimensions for V')
end

% JULIA ___________________________________________________________________
% add current path in Julia
a = jl.eval('LOAD_PATH');
if a{end}~='.'
    jleval push!(LOAD_PATH, pwd() * "\\external\\DpROM");
    jleval push!(LOAD_PATH, pwd() * "\\DpROM");
    jleval push!(LOAD_PATH, ".")
end
% load packages
jleval using TensorOperations
jleval using DpROM

fprintf(' dKdxi, assembling %d elements ...', nel)

% call the function in Julia to compute all the tensors
dKdxi = jl.call('stiffness_matrix_sensitivity', elements, ...
    nodes, conn, C, U, XGauss, WGauss, formulation);

fprintf(' %.2f s\n',toc(t0))




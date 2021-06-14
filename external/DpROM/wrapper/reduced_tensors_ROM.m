% reduced_tensors_ROM
%
% Synthax:
% tensors = reduced_tensors_ROM(myAssembly, elements, V)
%
% Description: this is a wrapper function for YetAnotherFEcode that
% computes the "standard" reduced order tensors, as described in the paper.
% INPUTS
%   - myAssembly: Assembly from YetAnotherFEcode.
%   - elements: table of the elements
%   - V: Reduced Order Basis (unconstrained)
% OUTPUT
%   tensors: a struct variable with the following fields*:
%       .Q2             n*n         reduced stiffness tensor
%    	.Q3             n*n*n       reduced stiffness tensor
%   	.Q4             n*n*n*n     reduced stiffness tensor
%      	.time           computational time
%   *being n=size(V,2)
%
% Additional notes:
%   - ALL the elements are assumed to have the same properties in terms
%     of MATERIAL and QUADRATURE rules.
%   - this function uses the red_stiff_tensors function implemented in the
%     Julia module "DpROM.jl". 
%   - as such, this function supports ONLY models meshed with the elements
%     supported by both YetAnotherFEcode AND the DpROM.jl
%   - List of currently supported elements: 
%     Q8, TET10, HEX20, WED15               (in YetAnotherFEcode)
%     Q8, TET10, HEX20, WED15, Q4, HEX8     (in DpROM.jl)
%
% Reference: J. Marconi, P. Tiso, D.E. Quadrelli & F. Braghin, "A higher 
% order parametric nonlinear reduced order model for imperfect structures 
% using Neumann expansion", Nonlinear Dynamics, 2021.
%
% Created: 14 May 2021
% Author: Jacopo Marconi, Politecnico di Milano

function tensors = reduced_tensors_ROM(myAssembly, elements, V)

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
if size(V,1)==nfree
    U = zeros(ndofs,size(V,2));
    U(freedofs,:) = V;
    V = U;
elseif size(V,1)~=ndofs
    error('Wrong dimensions for V')
end

% JULIA ___________________________________________________________________
% add current path in Julia
a = jl.eval('LOAD_PATH');
if a{end}~='.' % a{end}~='.'
    jleval push!(LOAD_PATH, pwd() * "\\external\\DpROM");
    jleval push!(LOAD_PATH, pwd() * "\\DpROM");
    jleval push!(LOAD_PATH, ".")
    fprintf(' Path added\n\n')
end
% load packages
jleval using TensorOperations
jleval using DpROM

% call the function once (for 1 element only) to precompile
elem1 = elements(1,:);
jl.call('red_stiff_tensors', elem1, ...
    nodes, conn, C, V, XGauss, WGauss);

disp(' REDUCED TENSORS (standard ~ Julia):')
fprintf(' Assembling %d elements ...', nel)

% call the function in Julia to compute all the tensors

%  a=jl.call('red_stiff_tensors', elements, ...
  %    nodes, conn, C, V, XGauss, WGauss);

%--------------------------------------------------- i added this
b=cell(1,10);
[b{1:10}]=jl.call('red_stiff_tensors', elements,nodes, conn, C, V, XGauss, WGauss);
    
%----------------------------------------------------

% unpack results __________________________________________________________
% Q2 = getfield(a,'1'); %#ok<*GFLD>
% Q3 = tensor(getfield(a,'2'));
% Q4 = tensor(getfield(a,'3'));
% time = double(getfield(a,'4'))/1000;

Q2=cell2mat(b(1,1));
Q3=tensor(cell2mat(b(1,2)));
Q4=tensor(cell2mat(b(1,3)));
time=double(cell2mat(b(1,4)))/1000;




fprintf(' %.2f s (%.2f s)\n',toc(t0),time)
fprintf(' SPEED: %.1f el/s\n',nel/time)
fprintf(' SIZEs: %d \n\n', size(V,2))

tensors.Q2 = Q2;
tensors.Q3 = Q3;
tensors.Q4 = Q4;
tensors.time = time;

% compute tensors for the tangent stiffness matrix (see tensors_KF.m)
tensors.Q3t = Q3 + permute(Q3, [1 3 2]); 
tensors.Q4t = Q4 + permute(Q4, [1 3 2 4]) + permute(Q4, [1 4 2 3]);

end



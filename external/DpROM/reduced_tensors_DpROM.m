% reduced_tensors_DpROM
%
% Synthax:
% tensors = reduced_tensors_DpROM(myAssembly, elements, V, U, FORMULATION, volume)
%
% Description: this is a wrapper function for YetAnotherFEcode that
% computes the parametric reduced order tensors for shape defect described
% in the paper.
% INPUTS
%   - myAssembly: Assembly from YetAnotherFEcode.
%   - elements: table of the elements
%   - V: Reduced Order Basis (unconstrained)
%   - U: tall matrix collecting by columns the (unconstrained) defect shapes
%   - formulation (optional): choose the Neumann formulation between "N1",
%     "N1T" and "N0".
%     The default value is N1.
%   - volume (optional): set to 1 to compute the additional tensors required
%     to approximate the volume integral over the defected volume, set to 0
%     if the integral can be computed over the nominal volume (e.g. 
%     isochoric defect).
%     The default value is 1.
% OUTPUT
%   tensors: a struct variable with the following fields*:
%       .Q2             n*n         reduced stiffness tensor
%    	.Q3             n*n*n       reduced stiffness tensor
%   	.Q4             n*n*n*n     reduced stiffness tensor
%    	.Q3d            n*n*d       reduced stiffness tensor
%       .Q4d            n*n*d*n     reduced stiffness tensor
%       .Q5d            n*n*n*d*n   reduced stiffness tensor
%      	.Q4dd           n*n*d*d     reduced stiffness tensor
%      	.Q5dd           n*n*d*n*d   reduced stiffness tensor
%     	.Q6dd           n*n*d*n*d*n reduced stiffness tensor
%      	.time           computational time
%      	.formulation    adopted formulation
%      	.volume         1/0
%   *being n=size(V,2) and d=size(U,2)
%
% Additional notes:
%   - ALL the elements are assumed to have the same properties in terms
%     of MATERIAL and QUADRATURE rules.
%   - this function uses the red_stiff_tensors_defects function 
%     implemented in the Julia module "DpROM.jl". 
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

function tensors = reduced_tensors_DpROM(myAssembly, elements, V, U, FORMULATION, volume)

t0=tic;

% Handle incomplete input
if nargin < 6
    volume = 1;
elseif nargin < 5
    volume = 1;
    FORMULATION = 'N1';
end

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
nd = size(U,2);                                 % number of defects

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
    u = zeros(ndofs,size(V,2));
    u(freedofs,:) = V;
    V = u;
elseif size(V,1)~=ndofs
    error('Wrong dimensions for V')
end
if size(U,1)==nfree
    u = zeros(ndofs,size(U,2));
    u(freedofs,:) = U;
    U = u;
elseif isscalar(U) && U == 0
    U = V(:,1)*0;
elseif size(U,1)~=ndofs
    error('Wrong dimensions for U')
end
FORMULATION = upper(FORMULATION);

% JULIA ___________________________________________________________________
% add current path in Julia
a = jl.eval('LOAD_PATH');
if a{end}~='.'
    jleval push!(LOAD_PATH, pwd() * "\\external\\DpROM");
    jleval push!(LOAD_PATH, ".")
    fprintf(' Path added\n\n')
end
% load packages
jleval using TensorOperations
jleval using DpROM

% call the function once (for 1 element only) to precompile
elem1 = elements(1,:);
jl.call('red_stiff_tensors_defects', FORMULATION, volume, elem1, ...
    nodes, conn, C, V, U, XGauss, WGauss);

disp([' REDUCED TENSORS (' FORMULATION ' ~ Julia):'])
fprintf(' Assembling %d elements ...', nel)

% call the function in Julia to compute all the tensors
a=jl.call('red_stiff_tensors_defects', FORMULATION, volume, elements, ...
    nodes, conn, C, V, U, XGauss, WGauss);

% unpack results __________________________________________________________
Q2n= getfield(a,'1'); %#ok<*GFLD>
a2 = getfield(a,'2');
a3 = getfield(a,'3');
a4 = getfield(a,'4');
a5 = getfield(a,'5');
a6 = getfield(a,'6');
a7 = getfield(a,'7');
a8 = getfield(a,'8');
a9 = getfield(a,'9');
if volume == 1
    ind = nd+1;
else
    ind = 1; % (the other tensors are zero)
end
for ii = 1 : ind
    Q3d{ii}  = tensor(a2{ii}); %#ok<*AGROW>
    Q4dd{ii} = tensor(a3{ii});
    Q3n{ii}  = tensor(a4{ii});
    Q4d{ii}  = tensor(a5{ii});
    Q5dd{ii} = tensor(a6{ii});
    Q4n{ii}  = tensor(a7{ii});
    Q5d{ii}  = tensor(a8{ii});
    Q6dd{ii} = tensor(a9{ii});
end

time = double(getfield(a,'10'))/1000;

fprintf(' %.2f s (%.2f s)\n',toc(t0),time)
fprintf(' SPEED: %.1f el/s\n',nel/time)
fprintf(' SIZEs: %d - %d \n\n', size(V,2), size(U,2))

tensors.Q2   = Q2n;
tensors.Q3   = Q3n;
tensors.Q4   = Q4n;
tensors.Q3d  = Q3d;
tensors.Q4d  = Q4d;
tensors.Q5d  = Q5d;
tensors.Q4dd = Q4dd;
tensors.Q5dd = Q5dd;
tensors.Q6dd = Q6dd;
tensors.time = time;
tensors.formulation = FORMULATION;
tensors.volume = volume;

end



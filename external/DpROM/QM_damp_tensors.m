
function CTEN = QM_damp_tensors(alpha, beta, myAssembly, Phi, Theta)

t0=tic;
fprintf(' DAMPING REDUCED TENSORS (Quadratic Manifold) ... ')

nel = myAssembly.Mesh.nElements;

n = size(Phi,1);
m = size(Phi,2);
C2 = zeros(m,m);
C3a = zeros(m,m,m);
C3b = zeros(m,m,m);
C4 = zeros(m,m,m,m);
u0 = zeros(n,1);
for ee = 1 : nel
    elem = myAssembly.Mesh.Elements(ee).Object;
    elem_dofs = elem.iDOFs;
    Ke = elem.tangent_stiffness_and_force(u0);
    Me = elem.mass_matrix;
    Ce = alpha*Me + beta*Ke;
    Phi_e = Phi(elem_dofs,:);
    Theta_e = Theta(elem_dofs,:,:);
    
    C2 = C2 + Phi_e'*Ce*Phi_e;
    C3a = C3a + einsum('iI,ij,jJK->IJK', Phi_e, Ce, Theta_e);
    C3b = C3b + einsum('iIJ,ij,jK->IJK', Theta_e, Ce, Phi_e);
    C4 = C4 + einsum('iIJ,ij,jKL->IJKL', Theta_e, Ce, Theta_e);
end

CTEN.C2 = C2;
CTEN.C3a = C3a;
CTEN.C3b = C3b;
CTEN.C4 = C4;

time = toc(t0);
CTEN.time = time;
fprintf(' %.2f s \n', time)
fprintf(' SPEED: %.1f el/s \n\n', nel/time)
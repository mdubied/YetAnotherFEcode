
function MTEN = QM_mass_tensors(myAssembly, Phi, Theta)

t0=tic;
fprintf(' MASS REDUCED TENSORS (Quadratic Manifold) ... ')

nel = myAssembly.Mesh.nElements;

m = size(Phi,2);
M2 = zeros(m,m);
M3a = zeros(m,m,m);
M3b = zeros(m,m,m);
M4 = zeros(m,m,m,m);
for ee = 1 : nel
    elem = myAssembly.Mesh.Elements(ee).Object;
    elem_dofs = elem.iDOFs;
    Me = full(elem.mass_matrix);
    Phi_e = Phi(elem_dofs,:);
    Theta_e = Theta(elem_dofs,:,:);
    
    M2 = M2 + Phi_e'*Me*Phi_e;
    M3a = M3a + einsum('iI,ij,jJK->IJK', Phi_e, Me, Theta_e);
    M3b = M3b + einsum('iIJ,ij,jK->IJK', Theta_e, Me, Phi_e);
    M4 = M4 + einsum('iIJ,ij,jKL->IJKL', Theta_e, Me, Theta_e);
end

MTEN.M2 = M2;
MTEN.M3a = M3a;
MTEN.M3b = M3b;
MTEN.M4 = M4;

time = toc(t0);
MTEN.time = time;
fprintf(' %.2f s \n', time)
fprintf(' SPEED: %.1f el/s \n\n', nel/time)

function [GREt,GREx,GREy,GREz] = error_GREm(u,ub,model,M)

if nargin<4
    M = 1;
    Mx = 1;
    My = 1;
    Mz = 1;
end

nt = size(u,2);
ndofs = model.Ndofs;
freedofs = model.dofs.free;
nf = length(freedofs);

if size(u,1)==nf
    U = zeros(ndofs,nt);
    U(freedofs,:) = u;
    u = U;
end
if size(ub,1)==nf
    U = zeros(ndofs,nt);
    U(freedofs,:) = ub;
    ub = U;
end
if size(M,1)==nf
    MM = zeros(ndofs);
    MM(freedofs,freedofs) = M;
    M = MM;
end

if size(M,1)==ndofs
    Mx = M(1:3:end,1:3:end);
    My = M(2:3:end,2:3:end);
    Mz = M(3:3:end,3:3:end);
end

nom_GREt = 0;
den_GREt = 0;
nom_GREx = 0;
den_GREx = 0;
nom_GREy = 0;
den_GREy = 0;
nom_GREz = 0;
den_GREz = 0;
for tt = 1:nt
    ut  =  u(:,tt);
    ubt = ub(:,tt);
    nom_GREt = nom_GREt + (ut-ubt).'*M*(ut-ubt);
    den_GREt = den_GREt + ut.'*M*ut;
    
    ux  =  u(1:3:end,tt);
    ubx = ub(1:3:end,tt);
    nom_GREx = nom_GREx + (ux-ubx).'*Mx*(ux-ubx);
    den_GREx = den_GREx + ux.'*Mx*ux;
    
    uy  =  u(2:3:end,tt);
    uby = ub(2:3:end,tt);
    nom_GREy = nom_GREy + (uy-uby).'*My*(uy-uby);
    den_GREy = den_GREy + uy.'*My*uy;
    
    uz  =  u(3:3:end,tt);
    ubz = ub(3:3:end,tt);
    nom_GREz = nom_GREz + (uz-ubz).'*Mz*(uz-ubz);
    den_GREz = den_GREz + uz.'*Mz*uz;
end
GREt = sqrt(nom_GREt/den_GREt)*100;
GREx = sqrt(nom_GREx/den_GREx)*100;
GREy = sqrt(nom_GREy/den_GREy)*100;
GREz = sqrt(nom_GREz/den_GREz)*100;
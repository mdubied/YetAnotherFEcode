
function [MDs,MDs_names,time] = MDs_computes(model,K0,VMs,NMDs,ng)

freedofs = model.dofs.free;

% check dimensions
if length(NMDs)==1
    NMDs = [1 NMDs];
end
if ~(size(VMs,1)==model.Ndofs)
    U = zeros(model.Ndofs,size(VMs,2));
    U(freedofs,:) = VMs;
    VMs = U;
end
if size(K0,1)==model.Ndofs
    K0 = K0(freedofs,freedofs);
end

VMs = VMs(:,NMDs(1):NMDs(end));

Nm = NMDs(end)-(NMDs(1)-1);
MDs = zeros(length(freedofs),Nm*(Nm+1)/2);
MDs_names = cell(size(MDs,2),1);

disp(' MODAL DERIVATIVES:')
t0 = tic;

count = 0;
for ii = 1:size(VMs,2)
    
    VM_i = VMs(:,ii);
    
    str = pad(sprintf(' dKdq_%d: ',ii),10);
    fprintf(str)
    dKdq_i = MDs_dKdq_assembly(model,VM_i,ng);
    t1 = tic;
    fprintf('\b, MDs computation ...')
    dKdq_i = dKdq_i(freedofs,freedofs);
    
    for jj = 1:size(VMs,2)
        if jj<ii
            continue
        end
        count = count + 1;
        MDs(:,count) = -K0\(dKdq_i*VMs(freedofs,jj));
        MDs_names{count} = sprintf('i=%d,j=%d',ii,jj);
    end
    fprintf(' %.2f s\n',toc(t1))
    
end
fprintf(' %d MDs computed in %.2f s \n\n',count,toc(t0))
time = toc(t0);
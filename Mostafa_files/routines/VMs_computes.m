%taken from Jac VM
function [VMs,f0,time] = VMs_computes(K,M,Nm,showeig)


% Compute eigenmodes
t0 = tic;
fprintf(' VIBRATION MODES ... ')
[Phi,D] = eigs(K,M,Nm,'SM');
[f0,ind] = sort(sqrt(diag(D))/2/pi);
VMs = Phi(:,ind);
time = toc(t0);
fprintf('%.2f\n\n',time)

% check on imaginary modes
if ~isreal(f0)
    warning(' Complex eigenfrequencies')
    disp(' ')
end

if showeig == 1
    disp(' Eigenfrequencies:')
    fprintf(' %.4f Hz \n',f0(1:Nm))
    disp(' ')
end


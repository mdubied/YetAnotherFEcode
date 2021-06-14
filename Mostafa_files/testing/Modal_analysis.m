%% modal analysis
%returns Ur which is the moadl amplitude and Ur_exp which is phi*eta
function [Ur Ur_exp w] =Modal_analysis(BeamAssembly,forceDOF,forceAmplitude,omega,tmax,h,Nmm)
nc=BeamAssembly.Mesh.EBC.constrainedDOFs(:,1)';
nf=BeamAssembly.Mesh.EBC.unconstrainedDOFs;
nt=length(nc)+length(nf);
M=BeamAssembly.DATA.M;
K=BeamAssembly.DATA.K;
C=BeamAssembly.DATA.C;

Mff=M(nf,nf);
Kff=K(nf,nf);
Cff=C(nf,nf);

Mfc=M(nf,nc);
Kfc=K(nf,nc);
Cfc=C(nf,nc);

Mcf=M(nc,nf);
Kcf=K(nc,nf);
Ccf=C(nc,nf);

Mcc=M(nc,nc);
Kcc=K(nc,nc);
Ccc=C(nc,nc);

% Nmm = 5;
[Phi,D] = eigs(Kff,Mff, Nmm, 'SM');
[f0,ind] = sort(sqrt(diag(D))/2/pi);
Phi = Phi(:,ind);
for ii = 1 : size(Phi, 2)
    
    Phi(:,ii) = Phi(:,ii)/max(sqrt(sum(Phi(:,ii).^2,2)));
    
    Phi(:,ii) = Phi(:,ii) / (Phi(:,ii)'*(Mff)*Phi(:,ii));
end
FF=zeros(nt,1);
FF(forceDOF)=forceAmplitude;

Ff=FF(nf);
Fc=FF(nc);



w = linspace(0,omega, 1000);

for ii = 1:length(w)
    wi = w(ii);
    
    % full response
%     Fc = -(Mfc*wi^2 + 1i*wi*Cfc + Kfc)*[ [zeros(28,1)]];
    U(:,ii) = (-Mff*wi^2 + 1i*wi*Cff + Kff)^-1 * Ff;
    
    % modal response
%     Fcr = Phi'*Fc;
Ffr = Phi'*Ff;
    Ur(:,ii) = (Phi'*(-Mff*wi^2 + 1i*wi*Cff + Kff)*Phi)^-1 * Ffr;
end
Ur_exp = Phi*Ur;
titolo = '';
figure
semilogy(w,abs(U(forceDOF,:)))
hold on
semilogy(w,abs(Ur_exp(forceDOF,:)),'r--')
title(sprintf('%d dofs, %d modes %s',size(Kff,1),Nmm,titolo))

 t=linspace(0,tmax,round(tmax/h));
 Ut=zeros(length(t),1);
 for ii=1:length(t)
 Ut(ii)=Ut(ii)+abs(U(forceDOF,:))*cos(omega*t(ii));
%  Urt=abs(Ur_exp)*cos(omega*t(ii));
 end
 figure
 plot(t,Ut(forceDOF,:))
%  hold on 
%  plot(t,Urt(forceDOF,:),'r--')
end
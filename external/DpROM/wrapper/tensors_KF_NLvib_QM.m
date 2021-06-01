
function [Kt,fi] = tensors_KF_NLvib_QM(Q3,Q4,Q3t,Q4t,Phi,Theta,M,D,fext,q,qd,qdd)
% NB: in NLvib, the nonlinear function must contain ONLY the nonlinear
% terms (i.e. WITHOUT the linear ones)
fel = ttsv(Q3,q,-1)  + ttsv(Q4,q,-1);
Kt = ttsv(Q3t, q, -2) + ttsv(Q4t,q,-2);

T = double(Phi + squeeze(ttt(Theta,tensor(q),3,1)));
t = double(squeeze(ttt(ttt(Theta,tensor(qd),3,1),tensor(qd),2,1)));
fmass = (T'*M*T - Phi'*M*Phi)*qdd + T'*M*t;
fdamp = (T'*D*T - Phi'*D*Phi)*qd;
g = T'*fext;

fi = fel + fmass + fdamp - g;
function [Kt,fi] = tensors_KF(Q2,Q3,Q4,Q3t,Q4t,q)

fi = Q2*q + ttsv(Q3,q,-1) + ttsv(Q4,q,-1);
Kt = Q2 + ttsv(Q3t, q, -2) + ttsv(Q4t,q,-2);
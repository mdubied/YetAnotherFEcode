function der = twoDof_forces_der(q,p)

%set input
q1 = q(1);
q2 = q(2);

k11 = p(1);
k13 = p(2);
k21 = p(3);
k23 = p(4);

% define derivatives
% dfdq (only nonlinear forces in f)
dfdq = reshape([k13.*q1.^2.*3.0+k23.*(q1-q2).^2.*3.0,k23.*(q1-q2).^2.*-...
    3.0,k23.*(q1-q2).^2.*-3.0,k23.*(q1-q2).^2.*3.0],[2,2]);

% dfdp (linear + nonlinear forces in f)
dfdp = reshape([q1,0.0,q1.^3,0.0,q1-q2,-q1+q2,(q1-q2).^3,...
       -(q1-q2).^3],[2,4]);

% dfdq2 (nonlinear forces, dfdq2 for linear forces is null)
dfdq2 = reshape([k13.*q1.*6.0+k23.*(q1.*2.0-q2.*2.0).*3.0,...
        k23.*(q1.*2.0-q2.*2.0).*-3.0,k23.*(q1.*2.0-q2.*2.0).*-...
        3.0,k23.*(q1.*2.0-q2.*2.0).*3.0,k23.*(q1.*2.0-q2.*2.0).*-...
        3.0,k23.*(q1.*2.0-q2.*2.0).*3.0,k23.*(q1.*2.0-q2.*2.0).*3.0,...
        k23.*(q1.*2.0-q2.*2.0).*-3.0],[2,2,2]);

% dfdp2 (linear + nonlinear forces)
dfdp2 = zeros(2,4,4);

% dfdpdq (linear + nonlinear forces)
dfdqdp = reshape([1.0,0.0,0.0,0.0,q1.^2.*3.0,0.0,0.0,0.0,1.0,-1.0,-1.0,...
         1.0,(q1-q2).^2.*3.0,(q1-q2).^2.*-3.0,(q1-q2).^2.*-3.0,...
         (q1-q2).^2.*3.0],[2,2,4]);
     
% generate output in struct array
der.dfdq = dfdq;
der.dfdp = dfdp;
der.dfdq2 = dfdq2;
der.dfdp2 = dfdp2;
der.dfdqdp = dfdqdp;

end
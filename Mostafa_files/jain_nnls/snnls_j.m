function [E, xi]=snnls_j(G,b,tau)
ne=size(G,2);
E=[];
Z=1:ne;
xi=zeros(ne,1);
disp('Now solving sparse non negative least squared problem ...')
while(norm(G*xi-b)>tau*norm(b))
    Mu=G.'*(b-G*xi);
   
    [~,index]= max(Mu);
    E=[E,index];
    Z=setdiff(Z,index);
%     while 1
%         zeta = zeros(ne,1);
%         zeta(E) = pinv(G(:,E))*b;
%         if all(zeta(E)>0)
%             xi = zeta;
%             break;
%         end
%         alpha = min( xi(E)./(xi(E) - zeta(E)));
%         xi = xi + alpha*(zeta-xi);
%         E = find(xi>0);
%         Z = setdiff(1:ne,E);
%     end
    zeta=nnls_j(G(:,E),b);
    xi(E)=zeta;
end
disp(norm(G*xi-b))
disp(tau*norm(b))
disp('Elements chosen : ')
disp(E')
end

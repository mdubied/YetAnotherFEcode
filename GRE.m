 function GRError =GRE(u,ut)   %u=[txdof]

tt=size(u,2);
tt1=size(ut,2);
if tt ~= tt1
    disp('wrong dimension, u and ut must be same time dimensions')
end
GRE_num=0;
GRE_den=0;
for ii=1:tt
    GRE_num=GRE_num +(u(:,ii)-ut(:,ii))'*(u(:,ii)-ut(:,ii));
    GRE_den=GRE_den +u(:,ii)'*u(:,ii);

end
GRError= 100*sqrt(GRE_num)/sqrt(GRE_den);
 end
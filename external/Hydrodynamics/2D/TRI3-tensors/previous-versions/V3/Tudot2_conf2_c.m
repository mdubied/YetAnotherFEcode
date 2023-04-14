function T = Tudot2_conf2_c(rho,c,vw1,vw2,inv,x01,x02,x03,x04,x05,x06)
T(1,:) = [0,0,0,0,0,0];
T(2,:) = [0,0,0,0,0,0];

T(3,:) = [0,0,(1/2).*(2.*vw1.^2.*(2.*vw2.*(x03+(-1) ...
.*x05)+(-2).*vw1.*(x04+(-1).*x06)).^2.*((x03+(-1).*x05).^2+(x04+( ...
-1).*x06).^2).*((4.*vw1.^2+4.*vw2.^2).*((x03+(-1).*x05).^2+(x04+( ...
-1).*x06).^2)).^(-3/2)+(-1/2).*(2.*vw2.*(x03+(-1).*x05)+(-2).* ...
vw1.*(x04+(-1).*x06)).^2.*((4.*vw1.^2+4.*vw2.^2).*((x03+(-1).*x05) ...
.^2+(x04+(-1).*x06).^2)).^(-1/2)+2.*vw1.*(2.*vw2.*(x03+(-1).*x05)+ ...
(-2).*vw1.*(x04+(-1).*x06)).*((4.*vw1.^2+4.*vw2.^2).*((x03+(-1).* ...
x05).^2+(x04+(-1).*x06).^2)).^(-1/2).*(x04+(-1).*x06)+(-1/6).*c.* ...
inv.*(4.*vw1.^2+4.*vw2.^2).*(1+(-1).*inv.^2.*(4.*vw1.^2+4.*vw2.^2) ...
.^(-1).*(2.*vw2.*(x03+(-1).*x05)+(-2).*vw1.*(x04+(-1).*x06)).^2.*( ...
(x03+(-1).*x05).^2+(x04+(-1).*x06).^2).^(-1)).^(-1/2).*(2.*inv.* ...
vw1.*(2.*vw2.*(x03+(-1).*x05)+(-2).*vw1.*(x04+(-1).*x06)).*((x03+( ...
-1).*x05).^2+(x04+(-1).*x06).^2).*((4.*vw1.^2+4.*vw2.^2).*((x03+( ...
-1).*x05).^2+(x04+(-1).*x06).^2)).^(-3/2)+inv.*((4.*vw1.^2+4.* ...
vw2.^2).*((x03+(-1).*x05).^2+(x04+(-1).*x06).^2)).^(-1/2).*(x04+( ...
-1).*x06)).*(x04+(-1).*x06).*acos(inv.*(2.*vw2.*(x03+(-1).*x05)+( ...
-2).*vw1.*(x04+(-1).*x06)).*((4.*vw1.^2+4.*vw2.^2).*((x03+(-1).* ...
x05).^2+(x04+(-1).*x06).^2)).^(-1/2))+(1/12).*c.*inv.*vw1.*(x04+( ...
-1).*x06).*(pi.^2+(-4).*acos(inv.*(2.*vw2.*(x03+(-1).*x05)+(-2).* ...
vw1.*(x04+(-1).*x06)).*((4.*vw1.^2+4.*vw2.^2).*((x03+(-1).*x05) ...
.^2+(x04+(-1).*x06).^2)).^(-1/2)).^2)),(1/2).*(2.*vw1.*vw2.*(2.* ...
vw2.*(x03+(-1).*x05)+(-2).*vw1.*(x04+(-1).*x06)).^2.*((x03+(-1).* ...
x05).^2+(x04+(-1).*x06).^2).*((4.*vw1.^2+4.*vw2.^2).*((x03+(-1).* ...
x05).^2+(x04+(-1).*x06).^2)).^(-3/2)+2.*vw1.*((-1).*x03+x05).*(2.* ...
vw2.*(x03+(-1).*x05)+(-2).*vw1.*(x04+(-1).*x06)).*((4.*vw1.^2+4.* ...
vw2.^2).*((x03+(-1).*x05).^2+(x04+(-1).*x06).^2)).^(-1/2)+(-1/6).* ...
c.*inv.*(4.*vw1.^2+4.*vw2.^2).*(1+(-1).*inv.^2.*(4.*vw1.^2+4.* ...
vw2.^2).^(-1).*(2.*vw2.*(x03+(-1).*x05)+(-2).*vw1.*(x04+(-1).*x06) ...
).^2.*((x03+(-1).*x05).^2+(x04+(-1).*x06).^2).^(-1)).^(-1/2).*(2.* ...
inv.*vw2.*(2.*vw2.*(x03+(-1).*x05)+(-2).*vw1.*(x04+(-1).*x06)).*(( ...
x03+(-1).*x05).^2+(x04+(-1).*x06).^2).*((4.*vw1.^2+4.*vw2.^2).*(( ...
x03+(-1).*x05).^2+(x04+(-1).*x06).^2)).^(-3/2)+inv.*((-1).*x03+ ...
x05).*((4.*vw1.^2+4.*vw2.^2).*((x03+(-1).*x05).^2+(x04+(-1).*x06) ...
.^2)).^(-1/2)).*(x04+(-1).*x06).*acos(inv.*(2.*vw2.*(x03+(-1).* ...
x05)+(-2).*vw1.*(x04+(-1).*x06)).*((4.*vw1.^2+4.*vw2.^2).*((x03+( ...
-1).*x05).^2+(x04+(-1).*x06).^2)).^(-1/2))+(1/12).*c.*inv.*vw2.*( ...
x04+(-1).*x06).*(pi.^2+(-4).*acos(inv.*(2.*vw2.*(x03+(-1).*x05)+( ...
-2).*vw1.*(x04+(-1).*x06)).*((4.*vw1.^2+4.*vw2.^2).*((x03+(-1).* ...
x05).^2+(x04+(-1).*x06).^2)).^(-1/2)).^2)),(1/2).*(2.*vw1.^2.*(2.* ...
vw2.*(x03+(-1).*x05)+(-2).*vw1.*(x04+(-1).*x06)).^2.*((x03+(-1).* ...
x05).^2+(x04+(-1).*x06).^2).*((4.*vw1.^2+4.*vw2.^2).*((x03+(-1).* ...
x05).^2+(x04+(-1).*x06).^2)).^(-3/2)+(-1/2).*(2.*vw2.*(x03+(-1).* ...
x05)+(-2).*vw1.*(x04+(-1).*x06)).^2.*((4.*vw1.^2+4.*vw2.^2).*(( ...
x03+(-1).*x05).^2+(x04+(-1).*x06).^2)).^(-1/2)+2.*vw1.*(2.*vw2.*( ...
x03+(-1).*x05)+(-2).*vw1.*(x04+(-1).*x06)).*((4.*vw1.^2+4.*vw2.^2) ...
.*((x03+(-1).*x05).^2+(x04+(-1).*x06).^2)).^(-1/2).*(x04+(-1).* ...
x06)+(-1/6).*c.*inv.*(4.*vw1.^2+4.*vw2.^2).*(1+(-1).*inv.^2.*(4.* ...
vw1.^2+4.*vw2.^2).^(-1).*(2.*vw2.*(x03+(-1).*x05)+(-2).*vw1.*(x04+ ...
(-1).*x06)).^2.*((x03+(-1).*x05).^2+(x04+(-1).*x06).^2).^(-1)).^( ...
-1/2).*(2.*inv.*vw1.*(2.*vw2.*(x03+(-1).*x05)+(-2).*vw1.*(x04+(-1) ...
.*x06)).*((x03+(-1).*x05).^2+(x04+(-1).*x06).^2).*((4.*vw1.^2+4.* ...
vw2.^2).*((x03+(-1).*x05).^2+(x04+(-1).*x06).^2)).^(-3/2)+inv.*(( ...
4.*vw1.^2+4.*vw2.^2).*((x03+(-1).*x05).^2+(x04+(-1).*x06).^2)).^( ...
-1/2).*(x04+(-1).*x06)).*(x04+(-1).*x06).*acos(inv.*(2.*vw2.*(x03+ ...
(-1).*x05)+(-2).*vw1.*(x04+(-1).*x06)).*((4.*vw1.^2+4.*vw2.^2).*(( ...
x03+(-1).*x05).^2+(x04+(-1).*x06).^2)).^(-1/2))+(1/12).*c.*inv.* ...
vw1.*(x04+(-1).*x06).*(pi.^2+(-4).*acos(inv.*(2.*vw2.*(x03+(-1).* ...
x05)+(-2).*vw1.*(x04+(-1).*x06)).*((4.*vw1.^2+4.*vw2.^2).*((x03+( ...
-1).*x05).^2+(x04+(-1).*x06).^2)).^(-1/2)).^2)),(1/2).*(2.*vw1.* ...
vw2.*(2.*vw2.*(x03+(-1).*x05)+(-2).*vw1.*(x04+(-1).*x06)).^2.*(( ...
x03+(-1).*x05).^2+(x04+(-1).*x06).^2).*((4.*vw1.^2+4.*vw2.^2).*(( ...
x03+(-1).*x05).^2+(x04+(-1).*x06).^2)).^(-3/2)+2.*vw1.*((-1).*x03+ ...
x05).*(2.*vw2.*(x03+(-1).*x05)+(-2).*vw1.*(x04+(-1).*x06)).*((4.* ...
vw1.^2+4.*vw2.^2).*((x03+(-1).*x05).^2+(x04+(-1).*x06).^2)).^( ...
-1/2)+(-1/6).*c.*inv.*(4.*vw1.^2+4.*vw2.^2).*(1+(-1).*inv.^2.*(4.* ...
vw1.^2+4.*vw2.^2).^(-1).*(2.*vw2.*(x03+(-1).*x05)+(-2).*vw1.*(x04+ ...
(-1).*x06)).^2.*((x03+(-1).*x05).^2+(x04+(-1).*x06).^2).^(-1)).^( ...
-1/2).*(2.*inv.*vw2.*(2.*vw2.*(x03+(-1).*x05)+(-2).*vw1.*(x04+(-1) ...
.*x06)).*((x03+(-1).*x05).^2+(x04+(-1).*x06).^2).*((4.*vw1.^2+4.* ...
vw2.^2).*((x03+(-1).*x05).^2+(x04+(-1).*x06).^2)).^(-3/2)+inv.*(( ...
-1).*x03+x05).*((4.*vw1.^2+4.*vw2.^2).*((x03+(-1).*x05).^2+(x04+( ...
-1).*x06).^2)).^(-1/2)).*(x04+(-1).*x06).*acos(inv.*(2.*vw2.*(x03+ ...
(-1).*x05)+(-2).*vw1.*(x04+(-1).*x06)).*((4.*vw1.^2+4.*vw2.^2).*(( ...
x03+(-1).*x05).^2+(x04+(-1).*x06).^2)).^(-1/2))+(1/12).*c.*inv.* ...
vw2.*(x04+(-1).*x06).*(pi.^2+(-4).*acos(inv.*(2.*vw2.*(x03+(-1).* ...
x05)+(-2).*vw1.*(x04+(-1).*x06)).*((4.*vw1.^2+4.*vw2.^2).*((x03+( ...
-1).*x05).^2+(x04+(-1).*x06).^2)).^(-1/2)).^2))];

T(4,:) = [0,0,(1/2).*(2.* ...
vw1.*vw2.*(2.*vw2.*(x03+(-1).*x05)+(-2).*vw1.*(x04+(-1).*x06)) ...
.^2.*((x03+(-1).*x05).^2+(x04+(-1).*x06).^2).*((4.*vw1.^2+4.* ...
vw2.^2).*((x03+(-1).*x05).^2+(x04+(-1).*x06).^2)).^(-3/2)+2.*vw2.* ...
(2.*vw2.*(x03+(-1).*x05)+(-2).*vw1.*(x04+(-1).*x06)).*((4.*vw1.^2+ ...
4.*vw2.^2).*((x03+(-1).*x05).^2+(x04+(-1).*x06).^2)).^(-1/2).*( ...
x04+(-1).*x06)+(1/6).*c.*inv.*(4.*vw1.^2+4.*vw2.^2).*(x03+(-1).* ...
x05).*(1+(-1).*inv.^2.*(4.*vw1.^2+4.*vw2.^2).^(-1).*(2.*vw2.*(x03+ ...
(-1).*x05)+(-2).*vw1.*(x04+(-1).*x06)).^2.*((x03+(-1).*x05).^2+( ...
x04+(-1).*x06).^2).^(-1)).^(-1/2).*(2.*inv.*vw1.*(2.*vw2.*(x03+( ...
-1).*x05)+(-2).*vw1.*(x04+(-1).*x06)).*((x03+(-1).*x05).^2+(x04+( ...
-1).*x06).^2).*((4.*vw1.^2+4.*vw2.^2).*((x03+(-1).*x05).^2+(x04+( ...
-1).*x06).^2)).^(-3/2)+inv.*((4.*vw1.^2+4.*vw2.^2).*((x03+(-1).* ...
x05).^2+(x04+(-1).*x06).^2)).^(-1/2).*(x04+(-1).*x06)).*acos(inv.* ...
(2.*vw2.*(x03+(-1).*x05)+(-2).*vw1.*(x04+(-1).*x06)).*((4.*vw1.^2+ ...
4.*vw2.^2).*((x03+(-1).*x05).^2+(x04+(-1).*x06).^2)).^(-1/2))+( ...
-1/12).*c.*inv.*vw1.*(x03+(-1).*x05).*(pi.^2+(-4).*acos(inv.*(2.* ...
vw2.*(x03+(-1).*x05)+(-2).*vw1.*(x04+(-1).*x06)).*((4.*vw1.^2+4.* ...
vw2.^2).*((x03+(-1).*x05).^2+(x04+(-1).*x06).^2)).^(-1/2)).^2)),( ...
1/2).*(2.*vw2.^2.*(2.*vw2.*(x03+(-1).*x05)+(-2).*vw1.*(x04+(-1).* ...
x06)).^2.*((x03+(-1).*x05).^2+(x04+(-1).*x06).^2).*((4.*vw1.^2+4.* ...
vw2.^2).*((x03+(-1).*x05).^2+(x04+(-1).*x06).^2)).^(-3/2)+2.*vw2.* ...
((-1).*x03+x05).*(2.*vw2.*(x03+(-1).*x05)+(-2).*vw1.*(x04+(-1).* ...
x06)).*((4.*vw1.^2+4.*vw2.^2).*((x03+(-1).*x05).^2+(x04+(-1).*x06) ...
.^2)).^(-1/2)+(-1/2).*(2.*vw2.*(x03+(-1).*x05)+(-2).*vw1.*(x04+( ...
-1).*x06)).^2.*((4.*vw1.^2+4.*vw2.^2).*((x03+(-1).*x05).^2+(x04+( ...
-1).*x06).^2)).^(-1/2)+(1/6).*c.*inv.*(4.*vw1.^2+4.*vw2.^2).*(x03+ ...
(-1).*x05).*(1+(-1).*inv.^2.*(4.*vw1.^2+4.*vw2.^2).^(-1).*(2.* ...
vw2.*(x03+(-1).*x05)+(-2).*vw1.*(x04+(-1).*x06)).^2.*((x03+(-1).* ...
x05).^2+(x04+(-1).*x06).^2).^(-1)).^(-1/2).*(2.*inv.*vw2.*(2.* ...
vw2.*(x03+(-1).*x05)+(-2).*vw1.*(x04+(-1).*x06)).*((x03+(-1).*x05) ...
.^2+(x04+(-1).*x06).^2).*((4.*vw1.^2+4.*vw2.^2).*((x03+(-1).*x05) ...
.^2+(x04+(-1).*x06).^2)).^(-3/2)+inv.*((-1).*x03+x05).*((4.* ...
vw1.^2+4.*vw2.^2).*((x03+(-1).*x05).^2+(x04+(-1).*x06).^2)).^( ...
-1/2)).*acos(inv.*(2.*vw2.*(x03+(-1).*x05)+(-2).*vw1.*(x04+(-1).* ...
x06)).*((4.*vw1.^2+4.*vw2.^2).*((x03+(-1).*x05).^2+(x04+(-1).*x06) ...
.^2)).^(-1/2))+(-1/12).*c.*inv.*vw2.*(x03+(-1).*x05).*(pi.^2+(-4) ...
.*acos(inv.*(2.*vw2.*(x03+(-1).*x05)+(-2).*vw1.*(x04+(-1).*x06)).* ...
((4.*vw1.^2+4.*vw2.^2).*((x03+(-1).*x05).^2+(x04+(-1).*x06).^2)) ...
.^(-1/2)).^2)),(1/2).*(2.*vw1.*vw2.*(2.*vw2.*(x03+(-1).*x05)+(-2) ...
.*vw1.*(x04+(-1).*x06)).^2.*((x03+(-1).*x05).^2+(x04+(-1).*x06) ...
.^2).*((4.*vw1.^2+4.*vw2.^2).*((x03+(-1).*x05).^2+(x04+(-1).*x06) ...
.^2)).^(-3/2)+2.*vw2.*(2.*vw2.*(x03+(-1).*x05)+(-2).*vw1.*(x04+( ...
-1).*x06)).*((4.*vw1.^2+4.*vw2.^2).*((x03+(-1).*x05).^2+(x04+(-1) ...
.*x06).^2)).^(-1/2).*(x04+(-1).*x06)+(1/6).*c.*inv.*(4.*vw1.^2+4.* ...
vw2.^2).*(x03+(-1).*x05).*(1+(-1).*inv.^2.*(4.*vw1.^2+4.*vw2.^2) ...
.^(-1).*(2.*vw2.*(x03+(-1).*x05)+(-2).*vw1.*(x04+(-1).*x06)).^2.*( ...
(x03+(-1).*x05).^2+(x04+(-1).*x06).^2).^(-1)).^(-1/2).*(2.*inv.* ...
vw1.*(2.*vw2.*(x03+(-1).*x05)+(-2).*vw1.*(x04+(-1).*x06)).*((x03+( ...
-1).*x05).^2+(x04+(-1).*x06).^2).*((4.*vw1.^2+4.*vw2.^2).*((x03+( ...
-1).*x05).^2+(x04+(-1).*x06).^2)).^(-3/2)+inv.*((4.*vw1.^2+4.* ...
vw2.^2).*((x03+(-1).*x05).^2+(x04+(-1).*x06).^2)).^(-1/2).*(x04+( ...
-1).*x06)).*acos(inv.*(2.*vw2.*(x03+(-1).*x05)+(-2).*vw1.*(x04+( ...
-1).*x06)).*((4.*vw1.^2+4.*vw2.^2).*((x03+(-1).*x05).^2+(x04+(-1) ...
.*x06).^2)).^(-1/2))+(-1/12).*c.*inv.*vw1.*(x03+(-1).*x05).*( ...
pi.^2+(-4).*acos(inv.*(2.*vw2.*(x03+(-1).*x05)+(-2).*vw1.*(x04+( ...
-1).*x06)).*((4.*vw1.^2+4.*vw2.^2).*((x03+(-1).*x05).^2+(x04+(-1) ...
.*x06).^2)).^(-1/2)).^2)),(1/2).*(2.*vw2.^2.*(2.*vw2.*(x03+(-1).* ...
x05)+(-2).*vw1.*(x04+(-1).*x06)).^2.*((x03+(-1).*x05).^2+(x04+(-1) ...
.*x06).^2).*((4.*vw1.^2+4.*vw2.^2).*((x03+(-1).*x05).^2+(x04+(-1) ...
.*x06).^2)).^(-3/2)+2.*vw2.*((-1).*x03+x05).*(2.*vw2.*(x03+(-1).* ...
x05)+(-2).*vw1.*(x04+(-1).*x06)).*((4.*vw1.^2+4.*vw2.^2).*((x03+( ...
-1).*x05).^2+(x04+(-1).*x06).^2)).^(-1/2)+(-1/2).*(2.*vw2.*(x03+( ...
-1).*x05)+(-2).*vw1.*(x04+(-1).*x06)).^2.*((4.*vw1.^2+4.*vw2.^2).* ...
((x03+(-1).*x05).^2+(x04+(-1).*x06).^2)).^(-1/2)+(1/6).*c.*inv.*( ...
4.*vw1.^2+4.*vw2.^2).*(x03+(-1).*x05).*(1+(-1).*inv.^2.*(4.* ...
vw1.^2+4.*vw2.^2).^(-1).*(2.*vw2.*(x03+(-1).*x05)+(-2).*vw1.*(x04+ ...
(-1).*x06)).^2.*((x03+(-1).*x05).^2+(x04+(-1).*x06).^2).^(-1)).^( ...
-1/2).*(2.*inv.*vw2.*(2.*vw2.*(x03+(-1).*x05)+(-2).*vw1.*(x04+(-1) ...
.*x06)).*((x03+(-1).*x05).^2+(x04+(-1).*x06).^2).*((4.*vw1.^2+4.* ...
vw2.^2).*((x03+(-1).*x05).^2+(x04+(-1).*x06).^2)).^(-3/2)+inv.*(( ...
-1).*x03+x05).*((4.*vw1.^2+4.*vw2.^2).*((x03+(-1).*x05).^2+(x04+( ...
-1).*x06).^2)).^(-1/2)).*acos(inv.*(2.*vw2.*(x03+(-1).*x05)+(-2).* ...
vw1.*(x04+(-1).*x06)).*((4.*vw1.^2+4.*vw2.^2).*((x03+(-1).*x05) ...
.^2+(x04+(-1).*x06).^2)).^(-1/2))+(-1/12).*c.*inv.*vw2.*(x03+(-1) ...
.*x05).*(pi.^2+(-4).*acos(inv.*(2.*vw2.*(x03+(-1).*x05)+(-2).* ...
vw1.*(x04+(-1).*x06)).*((4.*vw1.^2+4.*vw2.^2).*((x03+(-1).*x05) ...
.^2+(x04+(-1).*x06).^2)).^(-1/2)).^2))];

T(5,:) = T(3,:);
T(6,:) = T(4,:);
T = 0.5*rho*T;
end
function T = Tu2_conf1_c(rho,c,vw1,vw2,inv,x01,x02,x03,x04,x05,x06)
T(1,:)=[(1/2).*((-1).*vw1.*(4.*vw1.^2+4.*vw2.^2).*(x01+(-1).*x03).*(2.* ...
vw2.*(x01+(-1).*x03)+(-2).*vw1.*(x02+(-1).*x04)).^2.*((4.*vw1.^2+ ...
4.*vw2.^2).*((x01+(-1).*x03).^2+(x02+(-1).*x04).^2)).^(-3/2)+4.* ...
vw1.*vw2.*(2.*vw2.*(x01+(-1).*x03)+(-2).*vw1.*(x02+(-1).*x04)).*(( ...
4.*vw1.^2+4.*vw2.^2).*((x01+(-1).*x03).^2+(x02+(-1).*x04).^2)).^( ...
-1/2)+(-1/6).*c.*inv.*(4.*vw1.^2+4.*vw2.^2).*(1+(-1).*inv.^2.*(4.* ...
vw1.^2+4.*vw2.^2).^(-1).*(2.*vw2.*(x01+(-1).*x03)+(-2).*vw1.*(x02+ ...
(-1).*x04)).^2.*((x01+(-1).*x03).^2+(x02+(-1).*x04).^2).^(-1)).^( ...
-1/2).*((-1).*inv.*(4.*vw1.^2+4.*vw2.^2).*(x01+(-1).*x03).*(2.* ...
vw2.*(x01+(-1).*x03)+(-2).*vw1.*(x02+(-1).*x04)).*((4.*vw1.^2+4.* ...
vw2.^2).*((x01+(-1).*x03).^2+(x02+(-1).*x04).^2)).^(-3/2)+2.*inv.* ...
vw2.*((4.*vw1.^2+4.*vw2.^2).*((x01+(-1).*x03).^2+(x02+(-1).*x04) ...
.^2)).^(-1/2)).*(x02+(-1).*x04).*acos(inv.*(2.*vw2.*(x01+(-1).* ...
x03)+(-2).*vw1.*(x02+(-1).*x04)).*((4.*vw1.^2+4.*vw2.^2).*((x01+( ...
-1).*x03).^2+(x02+(-1).*x04).^2)).^(-1/2))),(1/2).*((-4).*vw1.^2.* ...
(2.*vw2.*(x01+(-1).*x03)+(-2).*vw1.*(x02+(-1).*x04)).*((4.*vw1.^2+ ...
4.*vw2.^2).*((x01+(-1).*x03).^2+(x02+(-1).*x04).^2)).^(-1/2)+(-1) ...
.*vw1.*(4.*vw1.^2+4.*vw2.^2).*(2.*vw2.*(x01+(-1).*x03)+(-2).*vw1.* ...
(x02+(-1).*x04)).^2.*((4.*vw1.^2+4.*vw2.^2).*((x01+(-1).*x03).^2+( ...
x02+(-1).*x04).^2)).^(-3/2).*(x02+(-1).*x04)+(-1/6).*c.*inv.*(4.* ...
vw1.^2+4.*vw2.^2).*(1+(-1).*inv.^2.*(4.*vw1.^2+4.*vw2.^2).^(-1).*( ...
2.*vw2.*(x01+(-1).*x03)+(-2).*vw1.*(x02+(-1).*x04)).^2.*((x01+(-1) ...
.*x03).^2+(x02+(-1).*x04).^2).^(-1)).^(-1/2).*((-2).*inv.*vw1.*(( ...
4.*vw1.^2+4.*vw2.^2).*((x01+(-1).*x03).^2+(x02+(-1).*x04).^2)).^( ...
-1/2)+(-1).*inv.*(4.*vw1.^2+4.*vw2.^2).*(2.*vw2.*(x01+(-1).*x03)+( ...
-2).*vw1.*(x02+(-1).*x04)).*((4.*vw1.^2+4.*vw2.^2).*((x01+(-1).* ...
x03).^2+(x02+(-1).*x04).^2)).^(-3/2).*(x02+(-1).*x04)).*(x02+(-1) ...
.*x04).*acos(inv.*(2.*vw2.*(x01+(-1).*x03)+(-2).*vw1.*(x02+(-1).* ...
x04)).*((4.*vw1.^2+4.*vw2.^2).*((x01+(-1).*x03).^2+(x02+(-1).*x04) ...
.^2)).^(-1/2))+(-1/48).*c.*inv.*(4.*vw1.^2+4.*vw2.^2).*(pi.^2+(-4) ...
.*acos(inv.*(2.*vw2.*(x01+(-1).*x03)+(-2).*vw1.*(x02+(-1).*x04)).* ...
((4.*vw1.^2+4.*vw2.^2).*((x01+(-1).*x03).^2+(x02+(-1).*x04).^2)) ...
.^(-1/2)).^2)),(1/2).*(vw1.*(4.*vw1.^2+4.*vw2.^2).*(x01+(-1).*x03) ...
.*(2.*vw2.*(x01+(-1).*x03)+(-2).*vw1.*(x02+(-1).*x04)).^2.*((4.* ...
vw1.^2+4.*vw2.^2).*((x01+(-1).*x03).^2+(x02+(-1).*x04).^2)).^( ...
-3/2)+(-4).*vw1.*vw2.*(2.*vw2.*(x01+(-1).*x03)+(-2).*vw1.*(x02+( ...
-1).*x04)).*((4.*vw1.^2+4.*vw2.^2).*((x01+(-1).*x03).^2+(x02+(-1) ...
.*x04).^2)).^(-1/2)+(-1/6).*c.*inv.*(4.*vw1.^2+4.*vw2.^2).*(1+(-1) ...
.*inv.^2.*(4.*vw1.^2+4.*vw2.^2).^(-1).*(2.*vw2.*(x01+(-1).*x03)+( ...
-2).*vw1.*(x02+(-1).*x04)).^2.*((x01+(-1).*x03).^2+(x02+(-1).*x04) ...
.^2).^(-1)).^(-1/2).*(inv.*(4.*vw1.^2+4.*vw2.^2).*(x01+(-1).*x03) ...
.*(2.*vw2.*(x01+(-1).*x03)+(-2).*vw1.*(x02+(-1).*x04)).*((4.* ...
vw1.^2+4.*vw2.^2).*((x01+(-1).*x03).^2+(x02+(-1).*x04).^2)).^( ...
-3/2)+(-2).*inv.*vw2.*((4.*vw1.^2+4.*vw2.^2).*((x01+(-1).*x03).^2+ ...
(x02+(-1).*x04).^2)).^(-1/2)).*(x02+(-1).*x04).*acos(inv.*(2.* ...
vw2.*(x01+(-1).*x03)+(-2).*vw1.*(x02+(-1).*x04)).*((4.*vw1.^2+4.* ...
vw2.^2).*((x01+(-1).*x03).^2+(x02+(-1).*x04).^2)).^(-1/2))),(1/2) ...
.*(4.*vw1.^2.*(2.*vw2.*(x01+(-1).*x03)+(-2).*vw1.*(x02+(-1).*x04)) ...
.*((4.*vw1.^2+4.*vw2.^2).*((x01+(-1).*x03).^2+(x02+(-1).*x04).^2)) ...
.^(-1/2)+vw1.*(4.*vw1.^2+4.*vw2.^2).*(2.*vw2.*(x01+(-1).*x03)+(-2) ...
.*vw1.*(x02+(-1).*x04)).^2.*((4.*vw1.^2+4.*vw2.^2).*((x01+(-1).* ...
x03).^2+(x02+(-1).*x04).^2)).^(-3/2).*(x02+(-1).*x04)+(-1/6).*c.* ...
inv.*(4.*vw1.^2+4.*vw2.^2).*(1+(-1).*inv.^2.*(4.*vw1.^2+4.*vw2.^2) ...
.^(-1).*(2.*vw2.*(x01+(-1).*x03)+(-2).*vw1.*(x02+(-1).*x04)).^2.*( ...
(x01+(-1).*x03).^2+(x02+(-1).*x04).^2).^(-1)).^(-1/2).*(2.*inv.* ...
vw1.*((4.*vw1.^2+4.*vw2.^2).*((x01+(-1).*x03).^2+(x02+(-1).*x04) ...
.^2)).^(-1/2)+inv.*(4.*vw1.^2+4.*vw2.^2).*(2.*vw2.*(x01+(-1).*x03)+(-2).*vw1.*(x02+(-1).*x04)).*((4.*vw1.^2+4.*vw2.^2).*((x01+(-1).* ...
x03).^2+(x02+(-1).*x04).^2)).^(-3/2).*(x02+(-1).*x04)).*(x02+(-1) ...
.*x04).*acos(inv.*(2.*vw2.*(x01+(-1).*x03)+(-2).*vw1.*(x02+(-1).* ...
x04)).*((4.*vw1.^2+4.*vw2.^2).*((x01+(-1).*x03).^2+(x02+(-1).*x04) ...
.^2)).^(-1/2))+(1/48).*c.*inv.*(4.*vw1.^2+4.*vw2.^2).*(pi.^2+(-4) ...
.*acos(inv.*(2.*vw2.*(x01+(-1).*x03)+(-2).*vw1.*(x02+(-1).*x04)).* ...
((4.*vw1.^2+4.*vw2.^2).*((x01+(-1).*x03).^2+(x02+(-1).*x04).^2)) ...
.^(-1/2)).^2)),0,0];

T(2,:)=[(1/2).*((-1).*vw2.*(4.*vw1.^2+4.*vw2.^2).*(x01+ ...
(-1).*x03).*(2.*vw2.*(x01+(-1).*x03)+(-2).*vw1.*(x02+(-1).*x04)) ...
.^2.*((4.*vw1.^2+4.*vw2.^2).*((x01+(-1).*x03).^2+(x02+(-1).*x04) ...
.^2)).^(-3/2)+4.*vw2.^2.*(2.*vw2.*(x01+(-1).*x03)+(-2).*vw1.*(x02+ ...
(-1).*x04)).*((4.*vw1.^2+4.*vw2.^2).*((x01+(-1).*x03).^2+(x02+(-1) ...
.*x04).^2)).^(-1/2)+(1/6).*c.*inv.*(4.*vw1.^2+4.*vw2.^2).*(x01+( ...
-1).*x03).*(1+(-1).*inv.^2.*(4.*vw1.^2+4.*vw2.^2).^(-1).*(2.*vw2.* ...
(x01+(-1).*x03)+(-2).*vw1.*(x02+(-1).*x04)).^2.*((x01+(-1).*x03) ...
.^2+(x02+(-1).*x04).^2).^(-1)).^(-1/2).*((-1).*inv.*(4.*vw1.^2+4.* ...
vw2.^2).*(x01+(-1).*x03).*(2.*vw2.*(x01+(-1).*x03)+(-2).*vw1.*( ...
x02+(-1).*x04)).*((4.*vw1.^2+4.*vw2.^2).*((x01+(-1).*x03).^2+(x02+ ...
(-1).*x04).^2)).^(-3/2)+2.*inv.*vw2.*((4.*vw1.^2+4.*vw2.^2).*(( ...
x01+(-1).*x03).^2+(x02+(-1).*x04).^2)).^(-1/2)).*acos(inv.*(2.* ...
vw2.*(x01+(-1).*x03)+(-2).*vw1.*(x02+(-1).*x04)).*((4.*vw1.^2+4.* ...
vw2.^2).*((x01+(-1).*x03).^2+(x02+(-1).*x04).^2)).^(-1/2))+(1/48) ...
.*c.*inv.*(4.*vw1.^2+4.*vw2.^2).*(pi.^2+(-4).*acos(inv.*(2.*vw2.*( ...
x01+(-1).*x03)+(-2).*vw1.*(x02+(-1).*x04)).*((4.*vw1.^2+4.*vw2.^2) ...
.*((x01+(-1).*x03).^2+(x02+(-1).*x04).^2)).^(-1/2)).^2)),(1/2).*(( ...
-4).*vw1.*vw2.*(2.*vw2.*(x01+(-1).*x03)+(-2).*vw1.*(x02+(-1).*x04) ...
).*((4.*vw1.^2+4.*vw2.^2).*((x01+(-1).*x03).^2+(x02+(-1).*x04).^2) ...
).^(-1/2)+(-1).*vw2.*(4.*vw1.^2+4.*vw2.^2).*(2.*vw2.*(x01+(-1).* ...
x03)+(-2).*vw1.*(x02+(-1).*x04)).^2.*((4.*vw1.^2+4.*vw2.^2).*(( ...
x01+(-1).*x03).^2+(x02+(-1).*x04).^2)).^(-3/2).*(x02+(-1).*x04)+( ...
1/6).*c.*inv.*(4.*vw1.^2+4.*vw2.^2).*(x01+(-1).*x03).*(1+(-1).* ...
inv.^2.*(4.*vw1.^2+4.*vw2.^2).^(-1).*(2.*vw2.*(x01+(-1).*x03)+(-2) ...
.*vw1.*(x02+(-1).*x04)).^2.*((x01+(-1).*x03).^2+(x02+(-1).*x04) ...
.^2).^(-1)).^(-1/2).*((-2).*inv.*vw1.*((4.*vw1.^2+4.*vw2.^2).*(( ...
x01+(-1).*x03).^2+(x02+(-1).*x04).^2)).^(-1/2)+(-1).*inv.*(4.* ...
vw1.^2+4.*vw2.^2).*(2.*vw2.*(x01+(-1).*x03)+(-2).*vw1.*(x02+(-1).* ...
x04)).*((4.*vw1.^2+4.*vw2.^2).*((x01+(-1).*x03).^2+(x02+(-1).*x04) ...
.^2)).^(-3/2).*(x02+(-1).*x04)).*acos(inv.*(2.*vw2.*(x01+(-1).* ...
x03)+(-2).*vw1.*(x02+(-1).*x04)).*((4.*vw1.^2+4.*vw2.^2).*((x01+( ...
-1).*x03).^2+(x02+(-1).*x04).^2)).^(-1/2))),(1/2).*(vw2.*(4.* ...
vw1.^2+4.*vw2.^2).*(x01+(-1).*x03).*(2.*vw2.*(x01+(-1).*x03)+(-2) ...
.*vw1.*(x02+(-1).*x04)).^2.*((4.*vw1.^2+4.*vw2.^2).*((x01+(-1).* ...
x03).^2+(x02+(-1).*x04).^2)).^(-3/2)+(-4).*vw2.^2.*(2.*vw2.*(x01+( ...
-1).*x03)+(-2).*vw1.*(x02+(-1).*x04)).*((4.*vw1.^2+4.*vw2.^2).*(( ...
x01+(-1).*x03).^2+(x02+(-1).*x04).^2)).^(-1/2)+(1/6).*c.*inv.*(4.* ...
vw1.^2+4.*vw2.^2).*(x01+(-1).*x03).*(1+(-1).*inv.^2.*(4.*vw1.^2+ ...
4.*vw2.^2).^(-1).*(2.*vw2.*(x01+(-1).*x03)+(-2).*vw1.*(x02+(-1).* ...
x04)).^2.*((x01+(-1).*x03).^2+(x02+(-1).*x04).^2).^(-1)).^(-1/2).* ...
(inv.*(4.*vw1.^2+4.*vw2.^2).*(x01+(-1).*x03).*(2.*vw2.*(x01+(-1).* ...
x03)+(-2).*vw1.*(x02+(-1).*x04)).*((4.*vw1.^2+4.*vw2.^2).*((x01+( ...
-1).*x03).^2+(x02+(-1).*x04).^2)).^(-3/2)+(-2).*inv.*vw2.*((4.* ...
vw1.^2+4.*vw2.^2).*((x01+(-1).*x03).^2+(x02+(-1).*x04).^2)).^( ...
-1/2)).*acos(inv.*(2.*vw2.*(x01+(-1).*x03)+(-2).*vw1.*(x02+(-1).* ...
x04)).*((4.*vw1.^2+4.*vw2.^2).*((x01+(-1).*x03).^2+(x02+(-1).*x04) ...
.^2)).^(-1/2))+(-1/48).*c.*inv.*(4.*vw1.^2+4.*vw2.^2).*(pi.^2+(-4) ...
.*acos(inv.*(2.*vw2.*(x01+(-1).*x03)+(-2).*vw1.*(x02+(-1).*x04)).* ...
((4.*vw1.^2+4.*vw2.^2).*((x01+(-1).*x03).^2+(x02+(-1).*x04).^2)) ...
.^(-1/2)).^2)),(1/2).*(4.*vw1.*vw2.*(2.*vw2.*(x01+(-1).*x03)+(-2) ...
.*vw1.*(x02+(-1).*x04)).*((4.*vw1.^2+4.*vw2.^2).*((x01+(-1).*x03) ...
.^2+(x02+(-1).*x04).^2)).^(-1/2)+vw2.*(4.*vw1.^2+4.*vw2.^2).*(2.* ...
vw2.*(x01+(-1).*x03)+(-2).*vw1.*(x02+(-1).*x04)).^2.*((4.*vw1.^2+ ...
4.*vw2.^2).*((x01+(-1).*x03).^2+(x02+(-1).*x04).^2)).^(-3/2).*( ...
x02+(-1).*x04)+(1/6).*c.*inv.*(4.*vw1.^2+4.*vw2.^2).*(x01+(-1).* ...
x03).*(1+(-1).*inv.^2.*(4.*vw1.^2+4.*vw2.^2).^(-1).*(2.*vw2.*(x01+ ...
(-1).*x03)+(-2).*vw1.*(x02+(-1).*x04)).^2.*((x01+(-1).*x03).^2+( ...
x02+(-1).*x04).^2).^(-1)).^(-1/2).*(2.*inv.*vw1.*((4.*vw1.^2+4.* ...
vw2.^2).*((x01+(-1).*x03).^2+(x02+(-1).*x04).^2)).^(-1/2)+inv.*( ...
4.*vw1.^2+4.*vw2.^2).*(2.*vw2.*(x01+(-1).*x03)+(-2).*vw1.*(x02+( ...
-1).*x04)).*((4.*vw1.^2+4.*vw2.^2).*((x01+(-1).*x03).^2+(x02+(-1) ...
.*x04).^2)).^(-3/2).*(x02+(-1).*x04)).*acos(inv.*(2.*vw2.*(x01+( ...
-1).*x03)+(-2).*vw1.*(x02+(-1).*x04)).*((4.*vw1.^2+4.*vw2.^2).*(( ...
x01+(-1).*x03).^2+(x02+(-1).*x04).^2)).^(-1/2))),0,0];

T(3,:) = T(1,:);
T(4,:) = T(2,:);
T(5,:) = [0,0,0,0,0,0];
T(6,:) = [0,0,0,0,0,0];

T = 0.5*rho*T;

end
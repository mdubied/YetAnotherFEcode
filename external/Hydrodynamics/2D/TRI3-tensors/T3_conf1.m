function T = T3_conf1(rho,x01,x02,x03,x04,x05,x06)
T(:,:,1) = [(-1).*rho.*(x01+(-1).*x03).^2.*((x01+(-1).*x03).^2+(x02+(-1).* ...
x04).^2).^(-5/2).*(x02+(-1).*x04).^2+(1/3).*rho.*((x01+(-1).*x03) ...
.^2+(x02+(-1).*x04).^2).^(-3/2).*(x02+(-1).*x04).^2,(2/3).*rho.*( ...
x01+(-1).*x03).*((x01+(-1).*x03).^2+(x02+(-1).*x04).^2).^(-3/2).*( ...
x02+(-1).*x04)+(-1).*rho.*(x01+(-1).*x03).*((x01+(-1).*x03).^2+( ...
x02+(-1).*x04).^2).^(-5/2).*(x02+(-1).*x04).^3,rho.*(x01+(-1).* ...
x03).^2.*((x01+(-1).*x03).^2+(x02+(-1).*x04).^2).^(-5/2).*(x02+( ...
-1).*x04).^2+(-1/3).*rho.*((x01+(-1).*x03).^2+(x02+(-1).*x04).^2) ...
.^(-3/2).*(x02+(-1).*x04).^2,(-2/3).*rho.*(x01+(-1).*x03).*((x01+( ...
-1).*x03).^2+(x02+(-1).*x04).^2).^(-3/2).*(x02+(-1).*x04)+rho.*( ...
x01+(-1).*x03).*((x01+(-1).*x03).^2+(x02+(-1).*x04).^2).^(-5/2).*( ...
x02+(-1).*x04).^3,0,0;(2/3).*rho.*(x01+(-1).*x03).*((x01+(-1).* ...
x03).^2+(x02+(-1).*x04).^2).^(-3/2).*(x02+(-1).*x04)+(-1).*rho.*( ...
x01+(-1).*x03).*((x01+(-1).*x03).^2+(x02+(-1).*x04).^2).^(-5/2).*( ...
x02+(-1).*x04).^3,(-2/3).*rho.*((x01+(-1).*x03).^2+(x02+(-1).*x04) ...
.^2).^(-1/2)+(5/3).*rho.*((x01+(-1).*x03).^2+(x02+(-1).*x04).^2) ...
.^(-3/2).*(x02+(-1).*x04).^2+(-1).*rho.*((x01+(-1).*x03).^2+(x02+( ...
-1).*x04).^2).^(-5/2).*(x02+(-1).*x04).^4,(-2/3).*rho.*(x01+(-1).* ...
x03).*((x01+(-1).*x03).^2+(x02+(-1).*x04).^2).^(-3/2).*(x02+(-1).* ...
x04)+rho.*(x01+(-1).*x03).*((x01+(-1).*x03).^2+(x02+(-1).*x04).^2) ...
.^(-5/2).*(x02+(-1).*x04).^3,(2/3).*rho.*((x01+(-1).*x03).^2+(x02+ ...
(-1).*x04).^2).^(-1/2)+(-5/3).*rho.*((x01+(-1).*x03).^2+(x02+(-1) ...
.*x04).^2).^(-3/2).*(x02+(-1).*x04).^2+rho.*((x01+(-1).*x03).^2+( ...
x02+(-1).*x04).^2).^(-5/2).*(x02+(-1).*x04).^4,0,0;rho.*(x01+(-1) ...
.*x03).^2.*((x01+(-1).*x03).^2+(x02+(-1).*x04).^2).^(-5/2).*(x02+( ...
-1).*x04).^2+(-1/3).*rho.*((x01+(-1).*x03).^2+(x02+(-1).*x04).^2) ...
.^(-3/2).*(x02+(-1).*x04).^2,(-2/3).*rho.*(x01+(-1).*x03).*((x01+( ...
-1).*x03).^2+(x02+(-1).*x04).^2).^(-3/2).*(x02+(-1).*x04)+rho.*( ...
x01+(-1).*x03).*((x01+(-1).*x03).^2+(x02+(-1).*x04).^2).^(-5/2).*( ...
x02+(-1).*x04).^3,(-1).*rho.*(x01+(-1).*x03).^2.*((x01+(-1).*x03) ...
.^2+(x02+(-1).*x04).^2).^(-5/2).*(x02+(-1).*x04).^2+(1/3).*rho.*(( ...
x01+(-1).*x03).^2+(x02+(-1).*x04).^2).^(-3/2).*(x02+(-1).*x04).^2, ...
(2/3).*rho.*(x01+(-1).*x03).*((x01+(-1).*x03).^2+(x02+(-1).*x04) ...
.^2).^(-3/2).*(x02+(-1).*x04)+(-1).*rho.*(x01+(-1).*x03).*((x01+( ...
-1).*x03).^2+(x02+(-1).*x04).^2).^(-5/2).*(x02+(-1).*x04).^3,0,0;( ...
-2/3).*rho.*(x01+(-1).*x03).*((x01+(-1).*x03).^2+(x02+(-1).*x04) ...
.^2).^(-3/2).*(x02+(-1).*x04)+rho.*(x01+(-1).*x03).*((x01+(-1).* ...
x03).^2+(x02+(-1).*x04).^2).^(-5/2).*(x02+(-1).*x04).^3,(2/3).* ...
rho.*((x01+(-1).*x03).^2+(x02+(-1).*x04).^2).^(-1/2)+(-5/3).*rho.* ...
((x01+(-1).*x03).^2+(x02+(-1).*x04).^2).^(-3/2).*(x02+(-1).*x04) ...
.^2+rho.*((x01+(-1).*x03).^2+(x02+(-1).*x04).^2).^(-5/2).*(x02+( ...
-1).*x04).^4,(2/3).*rho.*(x01+(-1).*x03).*((x01+(-1).*x03).^2+( ...
x02+(-1).*x04).^2).^(-3/2).*(x02+(-1).*x04)+(-1).*rho.*(x01+(-1).* ...
x03).*((x01+(-1).*x03).^2+(x02+(-1).*x04).^2).^(-5/2).*(x02+(-1).* ...
x04).^3,(-2/3).*rho.*((x01+(-1).*x03).^2+(x02+(-1).*x04).^2).^( ...
-1/2)+(5/3).*rho.*((x01+(-1).*x03).^2+(x02+(-1).*x04).^2).^(-3/2) ...
.*(x02+(-1).*x04).^2+(-1).*rho.*((x01+(-1).*x03).^2+(x02+(-1).* ...
x04).^2).^(-5/2).*(x02+(-1).*x04).^4,0,0;0,0,0,0,0,0;0,0,0,0,0,0];

T(:,:,2) = zeros(6,6);

T(:,:,3) = T(:,:,1);
T(:,:,4) = T(:,:,2);
T(:,:,5) = T(:,:,1);
T(:,:,6) = T(:,:,2);

end
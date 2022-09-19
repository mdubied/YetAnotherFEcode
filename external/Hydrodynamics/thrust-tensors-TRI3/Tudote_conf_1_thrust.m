function T = Tudote_conf_1_thrust(A,rho,vw1,vw2,r11,r12,r21,r22,x01,x02,x03,x04,x05,x06,ud1,ud2,ud3,ud4,ud5,ud6)  
    T = [(-1/6).*A.*rho.*(r11.*(ud1+x01)+(-1).*r11.*(ud3+x03)+r12.*(ud2+( ...
  -1).*ud4+x02+(-1).*x04)).*(abs(vw1).^2+abs(vw2).^2).*(abs(ud1+(-1) ...
  .*ud3+x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^( ...
  -1/2).*(1+(-16).*(4.*vw1.^2+4.*vw2.^2).^(-2).*(vw1.*(r11.*((-1).* ...
  ud1+ud3+(-1).*x01+x03)+r12.*((-1).*ud2+ud4+(-1).*x02+x04))+vw2.*( ...
  r21.*((-1).*ud1+ud3+(-1).*x01+x03)+r22.*((-1).*ud2+ud4+(-1).*x02+ ...
  x04))).^2.*(abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+ ...
  x02+(-1).*x04).^2).^(-1)).^(-1/2).*(2.*(4.*vw1.^2+4.*vw2.^2).^(-1) ...
  .*((-1).*r11.*((-1).*ud1+ud3+(-1).*x01+x03)+(-1).*r12.*((-1).*ud2+ ...
  ud4+(-1).*x02+x04)).*(abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+ ...
  (-1).*ud4+x02+(-1).*x04).^2).^(-1/2)+16.*vw1.*(4.*vw1.^2+4.* ...
  vw2.^2).^(-2).*(vw1.*(r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*(( ...
  -1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+ ...
  x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).*(abs(ud1+(-1).*ud3+x01+ ...
  (-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2)).*acos( ...
  4.*(4.*vw1.^2+4.*vw2.^2).^(-1).*(vw1.*(r11.*((-1).*ud1+ud3+(-1).* ...
  x01+x03)+r12.*((-1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1).* ...
  ud1+ud3+(-1).*x01+x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).*(abs( ...
  ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).*x04) ...
  .^2).^(-1/2))+(1/48).*A.*rho.*(r11.*(ud1+x01)+(-1).*r11.*(ud3+x03) ...
  +r12.*(ud2+(-1).*ud4+x02+(-1).*x04)).*abs(vw1).*(abs(ud1+(-1).* ...
  ud3+x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2) ...
  .*(pi.^2+(-4).*acos(4.*(4.*vw1.^2+4.*vw2.^2).^(-1).*(vw1.*(r11.*(( ...
  -1).*ud1+ud3+(-1).*x01+x03)+r12.*((-1).*ud2+ud4+(-1).*x02+x04))+ ...
  vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+x03)+r22.*((-1).*ud2+ud4+(-1) ...
  .*x02+x04))).*(abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+(-1).* ...
  ud4+x02+(-1).*x04).^2).^(-1/2)).^2).*Sign(vw1),(-1/6).*A.*rho.*( ...
  r11.*(ud1+x01)+(-1).*r11.*(ud3+x03)+r12.*(ud2+(-1).*ud4+x02+(-1).* ...
  x04)).*(abs(vw1).^2+abs(vw2).^2).*(abs(ud1+(-1).*ud3+x01+(-1).* ...
  x03).^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2).*(1+(-16).*( ...
  4.*vw1.^2+4.*vw2.^2).^(-2).*(vw1.*(r11.*((-1).*ud1+ud3+(-1).*x01+ ...
  x03)+r12.*((-1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1).*ud1+ ...
  ud3+(-1).*x01+x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).^2.*(abs( ...
  ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).*x04) ...
  .^2).^(-1)).^(-1/2).*(2.*(4.*vw1.^2+4.*vw2.^2).^(-1).*((-1).*r21.* ...
  ((-1).*ud1+ud3+(-1).*x01+x03)+(-1).*r22.*((-1).*ud2+ud4+(-1).*x02+ ...
  x04)).*(abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+ ...
  (-1).*x04).^2).^(-1/2)+16.*vw2.*(4.*vw1.^2+4.*vw2.^2).^(-2).*( ...
  vw1.*(r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*((-1).*ud2+ud4+(-1) ...
  .*x02+x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+x03)+r22.*((-1).* ...
  ud2+ud4+(-1).*x02+x04))).*(abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+ ...
  abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2)).*acos(4.*(4.*vw1.^2+ ...
  4.*vw2.^2).^(-1).*(vw1.*(r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*( ...
  (-1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+ ...
  x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).*(abs(ud1+(-1).*ud3+x01+ ...
  (-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2))+(1/48) ...
  .*A.*rho.*(r11.*(ud1+x01)+(-1).*r11.*(ud3+x03)+r12.*(ud2+(-1).* ...
  ud4+x02+(-1).*x04)).*abs(vw2).*(abs(ud1+(-1).*ud3+x01+(-1).*x03) ...
  .^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2).*(pi.^2+(-4).* ...
  acos(4.*(4.*vw1.^2+4.*vw2.^2).^(-1).*(vw1.*(r11.*((-1).*ud1+ud3+( ...
  -1).*x01+x03)+r12.*((-1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1) ...
  .*ud1+ud3+(-1).*x01+x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).*( ...
  abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).* ...
  x04).^2).^(-1/2)).^2).*Sign(vw2),(-1/6).*A.*rho.*(r11.*(ud1+x01)+( ...
  -1).*r11.*(ud3+x03)+r12.*(ud2+(-1).*ud4+x02+(-1).*x04)).*(abs(vw1) ...
  .^2+abs(vw2).^2).*(abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+( ...
  -1).*ud4+x02+(-1).*x04).^2).^(-1/2).*(1+(-16).*(4.*vw1.^2+4.* ...
  vw2.^2).^(-2).*(vw1.*(r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*(( ...
  -1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+ ...
  x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).^2.*(abs(ud1+(-1).*ud3+ ...
  x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1)).^( ...
  -1/2).*(2.*(4.*vw1.^2+4.*vw2.^2).^(-1).*((-1).*r11.*((-1).*ud1+ ...
  ud3+(-1).*x01+x03)+(-1).*r12.*((-1).*ud2+ud4+(-1).*x02+x04)).*( ...
  abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).* ...
  x04).^2).^(-1/2)+16.*vw1.*(4.*vw1.^2+4.*vw2.^2).^(-2).*(vw1.*( ...
  r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*((-1).*ud2+ud4+(-1).*x02+ ...
  x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+x03)+r22.*((-1).*ud2+ ...
  ud4+(-1).*x02+x04))).*(abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs( ...
  ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2)).*acos(4.*(4.*vw1.^2+4.* ...
  vw2.^2).^(-1).*(vw1.*(r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*(( ...
  -1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+ ...
  x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).*(abs(ud1+(-1).*ud3+x01+ ...
  (-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2))+(1/48) ...
  .*A.*rho.*(r11.*(ud1+x01)+(-1).*r11.*(ud3+x03)+r12.*(ud2+(-1).* ...
  ud4+x02+(-1).*x04)).*abs(vw1).*(abs(ud1+(-1).*ud3+x01+(-1).*x03) ...
  .^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2).*(pi.^2+(-4).* ...
  acos(4.*(4.*vw1.^2+4.*vw2.^2).^(-1).*(vw1.*(r11.*((-1).*ud1+ud3+( ...
  -1).*x01+x03)+r12.*((-1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1) ...
  .*ud1+ud3+(-1).*x01+x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).*( ...
  abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).* ...
  x04).^2).^(-1/2)).^2).*Sign(vw1),(-1/6).*A.*rho.*(r11.*(ud1+x01)+( ...
  -1).*r11.*(ud3+x03)+r12.*(ud2+(-1).*ud4+x02+(-1).*x04)).*(abs(vw1) ...
  .^2+abs(vw2).^2).*(abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+( ...
  -1).*ud4+x02+(-1).*x04).^2).^(-1/2).*(1+(-16).*(4.*vw1.^2+4.* ...
  vw2.^2).^(-2).*(vw1.*(r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*(( ...
  -1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+ ...
  x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).^2.*(abs(ud1+(-1).*ud3+ ...
  x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1)).^( ...
  -1/2).*(2.*(4.*vw1.^2+4.*vw2.^2).^(-1).*((-1).*r21.*((-1).*ud1+ ...
  ud3+(-1).*x01+x03)+(-1).*r22.*((-1).*ud2+ud4+(-1).*x02+x04)).*( ...
  abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).* ...
  x04).^2).^(-1/2)+16.*vw2.*(4.*vw1.^2+4.*vw2.^2).^(-2).*(vw1.*( ...
  r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*((-1).*ud2+ud4+(-1).*x02+ ...
  x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+x03)+r22.*((-1).*ud2+ ...
  ud4+(-1).*x02+x04))).*(abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs( ...
  ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2)).*acos(4.*(4.*vw1.^2+4.* ...
  vw2.^2).^(-1).*(vw1.*(r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*(( ...
  -1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+ ...
  x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).*(abs(ud1+(-1).*ud3+x01+ ...
  (-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2))+(1/48) ...
  .*A.*rho.*(r11.*(ud1+x01)+(-1).*r11.*(ud3+x03)+r12.*(ud2+(-1).* ...
  ud4+x02+(-1).*x04)).*abs(vw2).*(abs(ud1+(-1).*ud3+x01+(-1).*x03) ...
  .^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2).*(pi.^2+(-4).* ...
  acos(4.*(4.*vw1.^2+4.*vw2.^2).^(-1).*(vw1.*(r11.*((-1).*ud1+ud3+( ...
  -1).*x01+x03)+r12.*((-1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1) ...
  .*ud1+ud3+(-1).*x01+x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).*( ...
  abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).* ...
  x04).^2).^(-1/2)).^2).*Sign(vw2),0,0;(-1/6).*A.*rho.*(r21.*(ud1+ ...
  x01)+(-1).*r21.*(ud3+x03)+r22.*(ud2+(-1).*ud4+x02+(-1).*x04)).*( ...
  abs(vw1).^2+abs(vw2).^2).*(abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+ ...
  abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2).*(1+(-16).*(4.* ...
  vw1.^2+4.*vw2.^2).^(-2).*(vw1.*(r11.*((-1).*ud1+ud3+(-1).*x01+x03) ...
  +r12.*((-1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1).*ud1+ud3+( ...
  -1).*x01+x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).^2.*(abs(ud1+( ...
  -1).*ud3+x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^( ...
  -1)).^(-1/2).*(2.*(4.*vw1.^2+4.*vw2.^2).^(-1).*((-1).*r11.*((-1).* ...
  ud1+ud3+(-1).*x01+x03)+(-1).*r12.*((-1).*ud2+ud4+(-1).*x02+x04)).* ...
  (abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).* ...
  x04).^2).^(-1/2)+16.*vw1.*(4.*vw1.^2+4.*vw2.^2).^(-2).*(vw1.*( ...
  r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*((-1).*ud2+ud4+(-1).*x02+ ...
  x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+x03)+r22.*((-1).*ud2+ ...
  ud4+(-1).*x02+x04))).*(abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs( ...
  ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2)).*acos(4.*(4.*vw1.^2+4.* ...
  vw2.^2).^(-1).*(vw1.*(r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*(( ...
  -1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+ ...
  x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).*(abs(ud1+(-1).*ud3+x01+ ...
  (-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2))+(1/48) ...
  .*A.*rho.*(r21.*(ud1+x01)+(-1).*r21.*(ud3+x03)+r22.*(ud2+(-1).* ...
  ud4+x02+(-1).*x04)).*abs(vw1).*(abs(ud1+(-1).*ud3+x01+(-1).*x03) ...
  .^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2).*(pi.^2+(-4).* ...
  acos(4.*(4.*vw1.^2+4.*vw2.^2).^(-1).*(vw1.*(r11.*((-1).*ud1+ud3+( ...
  -1).*x01+x03)+r12.*((-1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1) ...
  .*ud1+ud3+(-1).*x01+x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).*( ...
  abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).* ...
  x04).^2).^(-1/2)).^2).*Sign(vw1),(-1/6).*A.*rho.*(r21.*(ud1+x01)+( ...
  -1).*r21.*(ud3+x03)+r22.*(ud2+(-1).*ud4+x02+(-1).*x04)).*(abs(vw1) ...
  .^2+abs(vw2).^2).*(abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+( ...
  -1).*ud4+x02+(-1).*x04).^2).^(-1/2).*(1+(-16).*(4.*vw1.^2+4.* ...
  vw2.^2).^(-2).*(vw1.*(r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*(( ...
  -1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+ ...
  x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).^2.*(abs(ud1+(-1).*ud3+ ...
  x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1)).^( ...
  -1/2).*(2.*(4.*vw1.^2+4.*vw2.^2).^(-1).*((-1).*r21.*((-1).*ud1+ ...
  ud3+(-1).*x01+x03)+(-1).*r22.*((-1).*ud2+ud4+(-1).*x02+x04)).*( ...
  abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).* ...
  x04).^2).^(-1/2)+16.*vw2.*(4.*vw1.^2+4.*vw2.^2).^(-2).*(vw1.*( ...
  r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*((-1).*ud2+ud4+(-1).*x02+ ...
  x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+x03)+r22.*((-1).*ud2+ ...
  ud4+(-1).*x02+x04))).*(abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs( ...
  ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2)).*acos(4.*(4.*vw1.^2+4.* ...
  vw2.^2).^(-1).*(vw1.*(r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*(( ...
  -1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+ ...
  x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).*(abs(ud1+(-1).*ud3+x01+ ...
  (-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2))+(1/48) ...
  .*A.*rho.*(r21.*(ud1+x01)+(-1).*r21.*(ud3+x03)+r22.*(ud2+(-1).* ...
  ud4+x02+(-1).*x04)).*abs(vw2).*(abs(ud1+(-1).*ud3+x01+(-1).*x03) ...
  .^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2).*(pi.^2+(-4).* ...
  acos(4.*(4.*vw1.^2+4.*vw2.^2).^(-1).*(vw1.*(r11.*((-1).*ud1+ud3+( ...
  -1).*x01+x03)+r12.*((-1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1) ...
  .*ud1+ud3+(-1).*x01+x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).*( ...
  abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).* ...
  x04).^2).^(-1/2)).^2).*Sign(vw2),(-1/6).*A.*rho.*(r21.*(ud1+x01)+( ...
  -1).*r21.*(ud3+x03)+r22.*(ud2+(-1).*ud4+x02+(-1).*x04)).*(abs(vw1) ...
  .^2+abs(vw2).^2).*(abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+( ...
  -1).*ud4+x02+(-1).*x04).^2).^(-1/2).*(1+(-16).*(4.*vw1.^2+4.* ...
  vw2.^2).^(-2).*(vw1.*(r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*(( ...
  -1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+ ...
  x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).^2.*(abs(ud1+(-1).*ud3+ ...
  x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1)).^( ...
  -1/2).*(2.*(4.*vw1.^2+4.*vw2.^2).^(-1).*((-1).*r11.*((-1).*ud1+ ...
  ud3+(-1).*x01+x03)+(-1).*r12.*((-1).*ud2+ud4+(-1).*x02+x04)).*( ...
  abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).* ...
  x04).^2).^(-1/2)+16.*vw1.*(4.*vw1.^2+4.*vw2.^2).^(-2).*(vw1.*( ...
  r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*((-1).*ud2+ud4+(-1).*x02+ ...
  x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+x03)+r22.*((-1).*ud2+ ...
  ud4+(-1).*x02+x04))).*(abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs( ...
  ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2)).*acos(4.*(4.*vw1.^2+4.* ...
  vw2.^2).^(-1).*(vw1.*(r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*(( ...
  -1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+ ...
  x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).*(abs(ud1+(-1).*ud3+x01+ ...
  (-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2))+(1/48) ...
  .*A.*rho.*(r21.*(ud1+x01)+(-1).*r21.*(ud3+x03)+r22.*(ud2+(-1).* ...
  ud4+x02+(-1).*x04)).*abs(vw1).*(abs(ud1+(-1).*ud3+x01+(-1).*x03) ...
  .^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2).*(pi.^2+(-4).* ...
  acos(4.*(4.*vw1.^2+4.*vw2.^2).^(-1).*(vw1.*(r11.*((-1).*ud1+ud3+( ...
  -1).*x01+x03)+r12.*((-1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1) ...
  .*ud1+ud3+(-1).*x01+x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).*( ...
  abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).* ...
  x04).^2).^(-1/2)).^2).*Sign(vw1),(-1/6).*A.*rho.*(r21.*(ud1+x01)+( ...
  -1).*r21.*(ud3+x03)+r22.*(ud2+(-1).*ud4+x02+(-1).*x04)).*(abs(vw1) ...
  .^2+abs(vw2).^2).*(abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+( ...
  -1).*ud4+x02+(-1).*x04).^2).^(-1/2).*(1+(-16).*(4.*vw1.^2+4.* ...
  vw2.^2).^(-2).*(vw1.*(r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*(( ...
  -1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+ ...
  x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).^2.*(abs(ud1+(-1).*ud3+ ...
  x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1)).^( ...
  -1/2).*(2.*(4.*vw1.^2+4.*vw2.^2).^(-1).*((-1).*r21.*((-1).*ud1+ ...
  ud3+(-1).*x01+x03)+(-1).*r22.*((-1).*ud2+ud4+(-1).*x02+x04)).*( ...
  abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).* ...
  x04).^2).^(-1/2)+16.*vw2.*(4.*vw1.^2+4.*vw2.^2).^(-2).*(vw1.*( ...
  r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*((-1).*ud2+ud4+(-1).*x02+ ...
  x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+x03)+r22.*((-1).*ud2+ ...
  ud4+(-1).*x02+x04))).*(abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs( ...
  ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2)).*acos(4.*(4.*vw1.^2+4.* ...
  vw2.^2).^(-1).*(vw1.*(r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*(( ...
  -1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+ ...
  x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).*(abs(ud1+(-1).*ud3+x01+ ...
  (-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2))+(1/48) ...
  .*A.*rho.*(r21.*(ud1+x01)+(-1).*r21.*(ud3+x03)+r22.*(ud2+(-1).* ...
  ud4+x02+(-1).*x04)).*abs(vw2).*(abs(ud1+(-1).*ud3+x01+(-1).*x03) ...
  .^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2).*(pi.^2+(-4).* ...
  acos(4.*(4.*vw1.^2+4.*vw2.^2).^(-1).*(vw1.*(r11.*((-1).*ud1+ud3+( ...
  -1).*x01+x03)+r12.*((-1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1) ...
  .*ud1+ud3+(-1).*x01+x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).*( ...
  abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).* ...
  x04).^2).^(-1/2)).^2).*Sign(vw2),0,0;(-1/6).*A.*rho.*(r11.*(ud1+ ...
  x01)+(-1).*r11.*(ud3+x03)+r12.*(ud2+(-1).*ud4+x02+(-1).*x04)).*( ...
  abs(vw1).^2+abs(vw2).^2).*(abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+ ...
  abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2).*(1+(-16).*(4.* ...
  vw1.^2+4.*vw2.^2).^(-2).*(vw1.*(r11.*((-1).*ud1+ud3+(-1).*x01+x03) ...
  +r12.*((-1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1).*ud1+ud3+( ...
  -1).*x01+x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).^2.*(abs(ud1+( ...
  -1).*ud3+x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^( ...
  -1)).^(-1/2).*(2.*(4.*vw1.^2+4.*vw2.^2).^(-1).*((-1).*r11.*((-1).* ...
  ud1+ud3+(-1).*x01+x03)+(-1).*r12.*((-1).*ud2+ud4+(-1).*x02+x04)).* ...
  (abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).* ...
  x04).^2).^(-1/2)+16.*vw1.*(4.*vw1.^2+4.*vw2.^2).^(-2).*(vw1.*( ...
  r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*((-1).*ud2+ud4+(-1).*x02+ ...
  x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+x03)+r22.*((-1).*ud2+ ...
  ud4+(-1).*x02+x04))).*(abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs( ...
  ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2)).*acos(4.*(4.*vw1.^2+4.* ...
  vw2.^2).^(-1).*(vw1.*(r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*(( ...
  -1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+ ...
  x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).*(abs(ud1+(-1).*ud3+x01+ ...
  (-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2))+(1/48) ...
  .*A.*rho.*(r11.*(ud1+x01)+(-1).*r11.*(ud3+x03)+r12.*(ud2+(-1).* ...
  ud4+x02+(-1).*x04)).*abs(vw1).*(abs(ud1+(-1).*ud3+x01+(-1).*x03) ...
  .^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2).*(pi.^2+(-4).* ...
  acos(4.*(4.*vw1.^2+4.*vw2.^2).^(-1).*(vw1.*(r11.*((-1).*ud1+ud3+( ...
  -1).*x01+x03)+r12.*((-1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1) ...
  .*ud1+ud3+(-1).*x01+x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).*( ...
  abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).* ...
  x04).^2).^(-1/2)).^2).*Sign(vw1),(-1/6).*A.*rho.*(r11.*(ud1+x01)+( ...
  -1).*r11.*(ud3+x03)+r12.*(ud2+(-1).*ud4+x02+(-1).*x04)).*(abs(vw1) ...
  .^2+abs(vw2).^2).*(abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+( ...
  -1).*ud4+x02+(-1).*x04).^2).^(-1/2).*(1+(-16).*(4.*vw1.^2+4.* ...
  vw2.^2).^(-2).*(vw1.*(r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*(( ...
  -1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+ ...
  x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).^2.*(abs(ud1+(-1).*ud3+ ...
  x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1)).^( ...
  -1/2).*(2.*(4.*vw1.^2+4.*vw2.^2).^(-1).*((-1).*r21.*((-1).*ud1+ ...
  ud3+(-1).*x01+x03)+(-1).*r22.*((-1).*ud2+ud4+(-1).*x02+x04)).*( ...
  abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).* ...
  x04).^2).^(-1/2)+16.*vw2.*(4.*vw1.^2+4.*vw2.^2).^(-2).*(vw1.*( ...
  r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*((-1).*ud2+ud4+(-1).*x02+ ...
  x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+x03)+r22.*((-1).*ud2+ ...
  ud4+(-1).*x02+x04))).*(abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs( ...
  ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2)).*acos(4.*(4.*vw1.^2+4.* ...
  vw2.^2).^(-1).*(vw1.*(r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*(( ...
  -1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+ ...
  x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).*(abs(ud1+(-1).*ud3+x01+ ...
  (-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2))+(1/48) ...
  .*A.*rho.*(r11.*(ud1+x01)+(-1).*r11.*(ud3+x03)+r12.*(ud2+(-1).* ...
  ud4+x02+(-1).*x04)).*abs(vw2).*(abs(ud1+(-1).*ud3+x01+(-1).*x03) ...
  .^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2).*(pi.^2+(-4).* ...
  acos(4.*(4.*vw1.^2+4.*vw2.^2).^(-1).*(vw1.*(r11.*((-1).*ud1+ud3+( ...
  -1).*x01+x03)+r12.*((-1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1) ...
  .*ud1+ud3+(-1).*x01+x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).*( ...
  abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).* ...
  x04).^2).^(-1/2)).^2).*Sign(vw2),(-1/6).*A.*rho.*(r11.*(ud1+x01)+( ...
  -1).*r11.*(ud3+x03)+r12.*(ud2+(-1).*ud4+x02+(-1).*x04)).*(abs(vw1) ...
  .^2+abs(vw2).^2).*(abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+( ...
  -1).*ud4+x02+(-1).*x04).^2).^(-1/2).*(1+(-16).*(4.*vw1.^2+4.* ...
  vw2.^2).^(-2).*(vw1.*(r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*(( ...
  -1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+ ...
  x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).^2.*(abs(ud1+(-1).*ud3+ ...
  x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1)).^( ...
  -1/2).*(2.*(4.*vw1.^2+4.*vw2.^2).^(-1).*((-1).*r11.*((-1).*ud1+ ...
  ud3+(-1).*x01+x03)+(-1).*r12.*((-1).*ud2+ud4+(-1).*x02+x04)).*( ...
  abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).* ...
  x04).^2).^(-1/2)+16.*vw1.*(4.*vw1.^2+4.*vw2.^2).^(-2).*(vw1.*( ...
  r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*((-1).*ud2+ud4+(-1).*x02+ ...
  x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+x03)+r22.*((-1).*ud2+ ...
  ud4+(-1).*x02+x04))).*(abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs( ...
  ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2)).*acos(4.*(4.*vw1.^2+4.* ...
  vw2.^2).^(-1).*(vw1.*(r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*(( ...
  -1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+ ...
  x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).*(abs(ud1+(-1).*ud3+x01+ ...
  (-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2))+(1/48) ...
  .*A.*rho.*(r11.*(ud1+x01)+(-1).*r11.*(ud3+x03)+r12.*(ud2+(-1).* ...
  ud4+x02+(-1).*x04)).*abs(vw1).*(abs(ud1+(-1).*ud3+x01+(-1).*x03) ...
  .^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2).*(pi.^2+(-4).* ...
  acos(4.*(4.*vw1.^2+4.*vw2.^2).^(-1).*(vw1.*(r11.*((-1).*ud1+ud3+( ...
  -1).*x01+x03)+r12.*((-1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1) ...
  .*ud1+ud3+(-1).*x01+x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).*( ...
  abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).* ...
  x04).^2).^(-1/2)).^2).*Sign(vw1),(-1/6).*A.*rho.*(r11.*(ud1+x01)+( ...
  -1).*r11.*(ud3+x03)+r12.*(ud2+(-1).*ud4+x02+(-1).*x04)).*(abs(vw1) ...
  .^2+abs(vw2).^2).*(abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+( ...
  -1).*ud4+x02+(-1).*x04).^2).^(-1/2).*(1+(-16).*(4.*vw1.^2+4.* ...
  vw2.^2).^(-2).*(vw1.*(r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*(( ...
  -1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+ ...
  x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).^2.*(abs(ud1+(-1).*ud3+ ...
  x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1)).^( ...
  -1/2).*(2.*(4.*vw1.^2+4.*vw2.^2).^(-1).*((-1).*r21.*((-1).*ud1+ ...
  ud3+(-1).*x01+x03)+(-1).*r22.*((-1).*ud2+ud4+(-1).*x02+x04)).*( ...
  abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).* ...
  x04).^2).^(-1/2)+16.*vw2.*(4.*vw1.^2+4.*vw2.^2).^(-2).*(vw1.*( ...
  r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*((-1).*ud2+ud4+(-1).*x02+ ...
  x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+x03)+r22.*((-1).*ud2+ ...
  ud4+(-1).*x02+x04))).*(abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs( ...
  ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2)).*acos(4.*(4.*vw1.^2+4.* ...
  vw2.^2).^(-1).*(vw1.*(r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*(( ...
  -1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+ ...
  x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).*(abs(ud1+(-1).*ud3+x01+ ...
  (-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2))+(1/48) ...
  .*A.*rho.*(r11.*(ud1+x01)+(-1).*r11.*(ud3+x03)+r12.*(ud2+(-1).* ...
  ud4+x02+(-1).*x04)).*abs(vw2).*(abs(ud1+(-1).*ud3+x01+(-1).*x03) ...
  .^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2).*(pi.^2+(-4).* ...
  acos(4.*(4.*vw1.^2+4.*vw2.^2).^(-1).*(vw1.*(r11.*((-1).*ud1+ud3+( ...
  -1).*x01+x03)+r12.*((-1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1) ...
  .*ud1+ud3+(-1).*x01+x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).*( ...
  abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).* ...
  x04).^2).^(-1/2)).^2).*Sign(vw2),0,0;(-1/6).*A.*rho.*(r21.*(ud1+ ...
  x01)+(-1).*r21.*(ud3+x03)+r22.*(ud2+(-1).*ud4+x02+(-1).*x04)).*( ...
  abs(vw1).^2+abs(vw2).^2).*(abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+ ...
  abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2).*(1+(-16).*(4.* ...
  vw1.^2+4.*vw2.^2).^(-2).*(vw1.*(r11.*((-1).*ud1+ud3+(-1).*x01+x03) ...
  +r12.*((-1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1).*ud1+ud3+( ...
  -1).*x01+x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).^2.*(abs(ud1+( ...
  -1).*ud3+x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^( ...
  -1)).^(-1/2).*(2.*(4.*vw1.^2+4.*vw2.^2).^(-1).*((-1).*r11.*((-1).* ...
  ud1+ud3+(-1).*x01+x03)+(-1).*r12.*((-1).*ud2+ud4+(-1).*x02+x04)).* ...
  (abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).* ...
  x04).^2).^(-1/2)+16.*vw1.*(4.*vw1.^2+4.*vw2.^2).^(-2).*(vw1.*( ...
  r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*((-1).*ud2+ud4+(-1).*x02+ ...
  x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+x03)+r22.*((-1).*ud2+ ...
  ud4+(-1).*x02+x04))).*(abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs( ...
  ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2)).*acos(4.*(4.*vw1.^2+4.* ...
  vw2.^2).^(-1).*(vw1.*(r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*(( ...
  -1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+ ...
  x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).*(abs(ud1+(-1).*ud3+x01+ ...
  (-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2))+(1/48) ...
  .*A.*rho.*(r21.*(ud1+x01)+(-1).*r21.*(ud3+x03)+r22.*(ud2+(-1).* ...
  ud4+x02+(-1).*x04)).*abs(vw1).*(abs(ud1+(-1).*ud3+x01+(-1).*x03) ...
  .^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2).*(pi.^2+(-4).* ...
  acos(4.*(4.*vw1.^2+4.*vw2.^2).^(-1).*(vw1.*(r11.*((-1).*ud1+ud3+( ...
  -1).*x01+x03)+r12.*((-1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1) ...
  .*ud1+ud3+(-1).*x01+x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).*( ...
  abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).* ...
  x04).^2).^(-1/2)).^2).*Sign(vw1),(-1/6).*A.*rho.*(r21.*(ud1+x01)+( ...
  -1).*r21.*(ud3+x03)+r22.*(ud2+(-1).*ud4+x02+(-1).*x04)).*(abs(vw1) ...
  .^2+abs(vw2).^2).*(abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+( ...
  -1).*ud4+x02+(-1).*x04).^2).^(-1/2).*(1+(-16).*(4.*vw1.^2+4.* ...
  vw2.^2).^(-2).*(vw1.*(r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*(( ...
  -1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+ ...
  x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).^2.*(abs(ud1+(-1).*ud3+ ...
  x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1)).^( ...
  -1/2).*(2.*(4.*vw1.^2+4.*vw2.^2).^(-1).*((-1).*r21.*((-1).*ud1+ ...
  ud3+(-1).*x01+x03)+(-1).*r22.*((-1).*ud2+ud4+(-1).*x02+x04)).*( ...
  abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).* ...
  x04).^2).^(-1/2)+16.*vw2.*(4.*vw1.^2+4.*vw2.^2).^(-2).*(vw1.*( ...
  r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*((-1).*ud2+ud4+(-1).*x02+ ...
  x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+x03)+r22.*((-1).*ud2+ ...
  ud4+(-1).*x02+x04))).*(abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs( ...
  ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2)).*acos(4.*(4.*vw1.^2+4.* ...
  vw2.^2).^(-1).*(vw1.*(r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*(( ...
  -1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+ ...
  x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).*(abs(ud1+(-1).*ud3+x01+ ...
  (-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2))+(1/48) ...
  .*A.*rho.*(r21.*(ud1+x01)+(-1).*r21.*(ud3+x03)+r22.*(ud2+(-1).* ...
  ud4+x02+(-1).*x04)).*abs(vw2).*(abs(ud1+(-1).*ud3+x01+(-1).*x03) ...
  .^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2).*(pi.^2+(-4).* ...
  acos(4.*(4.*vw1.^2+4.*vw2.^2).^(-1).*(vw1.*(r11.*((-1).*ud1+ud3+( ...
  -1).*x01+x03)+r12.*((-1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1) ...
  .*ud1+ud3+(-1).*x01+x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).*( ...
  abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).* ...
  x04).^2).^(-1/2)).^2).*Sign(vw2),(-1/6).*A.*rho.*(r21.*(ud1+x01)+( ...
  -1).*r21.*(ud3+x03)+r22.*(ud2+(-1).*ud4+x02+(-1).*x04)).*(abs(vw1) ...
  .^2+abs(vw2).^2).*(abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+( ...
  -1).*ud4+x02+(-1).*x04).^2).^(-1/2).*(1+(-16).*(4.*vw1.^2+4.* ...
  vw2.^2).^(-2).*(vw1.*(r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*(( ...
  -1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+ ...
  x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).^2.*(abs(ud1+(-1).*ud3+ ...
  x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1)).^( ...
  -1/2).*(2.*(4.*vw1.^2+4.*vw2.^2).^(-1).*((-1).*r11.*((-1).*ud1+ ...
  ud3+(-1).*x01+x03)+(-1).*r12.*((-1).*ud2+ud4+(-1).*x02+x04)).*( ...
  abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).* ...
  x04).^2).^(-1/2)+16.*vw1.*(4.*vw1.^2+4.*vw2.^2).^(-2).*(vw1.*( ...
  r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*((-1).*ud2+ud4+(-1).*x02+ ...
  x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+x03)+r22.*((-1).*ud2+ ...
  ud4+(-1).*x02+x04))).*(abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs( ...
  ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2)).*acos(4.*(4.*vw1.^2+4.* ...
  vw2.^2).^(-1).*(vw1.*(r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*(( ...
  -1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+ ...
  x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).*(abs(ud1+(-1).*ud3+x01+ ...
  (-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2))+(1/48) ...
  .*A.*rho.*(r21.*(ud1+x01)+(-1).*r21.*(ud3+x03)+r22.*(ud2+(-1).* ...
  ud4+x02+(-1).*x04)).*abs(vw1).*(abs(ud1+(-1).*ud3+x01+(-1).*x03) ...
  .^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2).*(pi.^2+(-4).* ...
  acos(4.*(4.*vw1.^2+4.*vw2.^2).^(-1).*(vw1.*(r11.*((-1).*ud1+ud3+( ...
  -1).*x01+x03)+r12.*((-1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1) ...
  .*ud1+ud3+(-1).*x01+x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).*( ...
  abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).* ...
  x04).^2).^(-1/2)).^2).*Sign(vw1),(-1/6).*A.*rho.*(r21.*(ud1+x01)+( ...
  -1).*r21.*(ud3+x03)+r22.*(ud2+(-1).*ud4+x02+(-1).*x04)).*(abs(vw1) ...
  .^2+abs(vw2).^2).*(abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+( ...
  -1).*ud4+x02+(-1).*x04).^2).^(-1/2).*(1+(-16).*(4.*vw1.^2+4.* ...
  vw2.^2).^(-2).*(vw1.*(r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*(( ...
  -1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+ ...
  x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).^2.*(abs(ud1+(-1).*ud3+ ...
  x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1)).^( ...
  -1/2).*(2.*(4.*vw1.^2+4.*vw2.^2).^(-1).*((-1).*r21.*((-1).*ud1+ ...
  ud3+(-1).*x01+x03)+(-1).*r22.*((-1).*ud2+ud4+(-1).*x02+x04)).*( ...
  abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).* ...
  x04).^2).^(-1/2)+16.*vw2.*(4.*vw1.^2+4.*vw2.^2).^(-2).*(vw1.*( ...
  r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*((-1).*ud2+ud4+(-1).*x02+ ...
  x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+x03)+r22.*((-1).*ud2+ ...
  ud4+(-1).*x02+x04))).*(abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs( ...
  ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2)).*acos(4.*(4.*vw1.^2+4.* ...
  vw2.^2).^(-1).*(vw1.*(r11.*((-1).*ud1+ud3+(-1).*x01+x03)+r12.*(( ...
  -1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1).*ud1+ud3+(-1).*x01+ ...
  x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).*(abs(ud1+(-1).*ud3+x01+ ...
  (-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2))+(1/48) ...
  .*A.*rho.*(r21.*(ud1+x01)+(-1).*r21.*(ud3+x03)+r22.*(ud2+(-1).* ...
  ud4+x02+(-1).*x04)).*abs(vw2).*(abs(ud1+(-1).*ud3+x01+(-1).*x03) ...
  .^2+abs(ud2+(-1).*ud4+x02+(-1).*x04).^2).^(-1/2).*(pi.^2+(-4).* ...
  acos(4.*(4.*vw1.^2+4.*vw2.^2).^(-1).*(vw1.*(r11.*((-1).*ud1+ud3+( ...
  -1).*x01+x03)+r12.*((-1).*ud2+ud4+(-1).*x02+x04))+vw2.*(r21.*((-1) ...
  .*ud1+ud3+(-1).*x01+x03)+r22.*((-1).*ud2+ud4+(-1).*x02+x04))).*( ...
  abs(ud1+(-1).*ud3+x01+(-1).*x03).^2+abs(ud2+(-1).*ud4+x02+(-1).* ...
  x04).^2).^(-1/2)).^2).*Sign(vw2),0,0;0,0,0,0,0,0;0,0,0,0,0,0];
end

function T = Tudote_conf_3_drag(A,rho,vw1,vw2,r11,r12,r21,r22,x01,x02,x03,x04,x05,x06,ud1,ud2,ud3,ud4,ud5,ud6)
    T = [384.*A.*rho.*vw1.^2.*(4.*vw1.^2+4.*vw2.^2).^(-4).*(vw1.*(r11.*( ...
  ud1+(-1).*ud5+x01+(-1).*x05)+r12.*(ud2+(-1).*ud6+x02+(-1).*x06))+ ...
  vw2.*(r21.*(ud1+(-1).*ud5+x01+(-1).*x05)+r22.*(ud2+(-1).*ud6+x02+( ...
  -1).*x06))).^2.*(abs(vw1).^2+abs(vw2).^2).*(abs(ud1+(-1).*ud5+x01+ ...
  (-1).*x05).^2+abs(ud2+(-1).*ud6+x02+(-1).*x06).^2).^(-1)+(-16).* ...
  A.*rho.*(4.*vw1.^2+4.*vw2.^2).^(-3).*(vw1.*(r11.*(ud1+(-1).*ud5+ ...
  x01+(-1).*x05)+r12.*(ud2+(-1).*ud6+x02+(-1).*x06))+vw2.*(r21.*( ...
  ud1+(-1).*ud5+x01+(-1).*x05)+r22.*(ud2+(-1).*ud6+x02+(-1).*x06))) ...
  .^2.*(abs(vw1).^2+abs(vw2).^2).*(abs(ud1+(-1).*ud5+x01+(-1).*x05) ...
  .^2+abs(ud2+(-1).*ud6+x02+(-1).*x06).^2).^(-1)+32.*A.*rho.*vw1.*( ...
  4.*vw1.^2+4.*vw2.^2).^(-3).*(vw1.*(r11.*(ud1+(-1).*ud5+x01+(-1).* ...
  x05)+r12.*(ud2+(-1).*ud6+x02+(-1).*x06))+vw2.*(r21.*(ud1+(-1).* ...
  ud5+x01+(-1).*x05)+r22.*(ud2+(-1).*ud6+x02+(-1).*x06))).*((-1).* ...
  r11.*(ud1+(-1).*ud5+x01+(-1).*x05)+(-1).*r12.*(ud2+(-1).*ud6+x02+( ...
  -1).*x06)).*(abs(vw1).^2+abs(vw2).^2).*(abs(ud1+(-1).*ud5+x01+(-1) ...
  .*x05).^2+abs(ud2+(-1).*ud6+x02+(-1).*x06).^2).^(-1)+(-32).*A.* ...
  rho.*vw1.*(4.*vw1.^2+4.*vw2.^2).^(-3).*(vw1.*(r11.*(ud1+(-1).*ud5+ ...
  x01+(-1).*x05)+r12.*(ud2+(-1).*ud6+x02+(-1).*x06))+vw2.*(r21.*( ...
  ud1+(-1).*ud5+x01+(-1).*x05)+r22.*(ud2+(-1).*ud6+x02+(-1).*x06))) ...
  .^2.*abs(vw1).*(abs(ud1+(-1).*ud5+x01+(-1).*x05).^2+abs(ud2+(-1).* ...
  ud6+x02+(-1).*x06).^2).^(-1).*sign(vw1),384.*A.*rho.*vw1.*vw2.*( ...
  4.*vw1.^2+4.*vw2.^2).^(-4).*(vw1.*(r11.*(ud1+(-1).*ud5+x01+(-1).* ...
  x05)+r12.*(ud2+(-1).*ud6+x02+(-1).*x06))+vw2.*(r21.*(ud1+(-1).* ...
  ud5+x01+(-1).*x05)+r22.*(ud2+(-1).*ud6+x02+(-1).*x06))).^2.*(abs( ...
  vw1).^2+abs(vw2).^2).*(abs(ud1+(-1).*ud5+x01+(-1).*x05).^2+abs( ...
  ud2+(-1).*ud6+x02+(-1).*x06).^2).^(-1)+32.*A.*rho.*vw1.*(4.* ...
  vw1.^2+4.*vw2.^2).^(-3).*(vw1.*(r11.*(ud1+(-1).*ud5+x01+(-1).*x05) ...
  +r12.*(ud2+(-1).*ud6+x02+(-1).*x06))+vw2.*(r21.*(ud1+(-1).*ud5+ ...
  x01+(-1).*x05)+r22.*(ud2+(-1).*ud6+x02+(-1).*x06))).*((-1).*r21.*( ...
  ud1+(-1).*ud5+x01+(-1).*x05)+(-1).*r22.*(ud2+(-1).*ud6+x02+(-1).* ...
  x06)).*(abs(vw1).^2+abs(vw2).^2).*(abs(ud1+(-1).*ud5+x01+(-1).* ...
  x05).^2+abs(ud2+(-1).*ud6+x02+(-1).*x06).^2).^(-1)+(-32).*A.*rho.* ...
  vw1.*(4.*vw1.^2+4.*vw2.^2).^(-3).*(vw1.*(r11.*(ud1+(-1).*ud5+x01+( ...
  -1).*x05)+r12.*(ud2+(-1).*ud6+x02+(-1).*x06))+vw2.*(r21.*(ud1+(-1) ...
  .*ud5+x01+(-1).*x05)+r22.*(ud2+(-1).*ud6+x02+(-1).*x06))).^2.*abs( ...
  vw2).*(abs(ud1+(-1).*ud5+x01+(-1).*x05).^2+abs(ud2+(-1).*ud6+x02+( ...
  -1).*x06).^2).^(-1).*sign(vw2),0,0,384.*A.*rho.*vw1.^2.*(4.* ...
  vw1.^2+4.*vw2.^2).^(-4).*(vw1.*(r11.*(ud1+(-1).*ud5+x01+(-1).*x05) ...
  +r12.*(ud2+(-1).*ud6+x02+(-1).*x06))+vw2.*(r21.*(ud1+(-1).*ud5+ ...
  x01+(-1).*x05)+r22.*(ud2+(-1).*ud6+x02+(-1).*x06))).^2.*(abs(vw1) ...
  .^2+abs(vw2).^2).*(abs(ud1+(-1).*ud5+x01+(-1).*x05).^2+abs(ud2+( ...
  -1).*ud6+x02+(-1).*x06).^2).^(-1)+(-16).*A.*rho.*(4.*vw1.^2+4.* ...
  vw2.^2).^(-3).*(vw1.*(r11.*(ud1+(-1).*ud5+x01+(-1).*x05)+r12.*( ...
  ud2+(-1).*ud6+x02+(-1).*x06))+vw2.*(r21.*(ud1+(-1).*ud5+x01+(-1).* ...
  x05)+r22.*(ud2+(-1).*ud6+x02+(-1).*x06))).^2.*(abs(vw1).^2+abs( ...
  vw2).^2).*(abs(ud1+(-1).*ud5+x01+(-1).*x05).^2+abs(ud2+(-1).*ud6+ ...
  x02+(-1).*x06).^2).^(-1)+32.*A.*rho.*vw1.*(4.*vw1.^2+4.*vw2.^2).^( ...
  -3).*(vw1.*(r11.*(ud1+(-1).*ud5+x01+(-1).*x05)+r12.*(ud2+(-1).* ...
  ud6+x02+(-1).*x06))+vw2.*(r21.*(ud1+(-1).*ud5+x01+(-1).*x05)+r22.* ...
  (ud2+(-1).*ud6+x02+(-1).*x06))).*((-1).*r11.*(ud1+(-1).*ud5+x01+( ...
  -1).*x05)+(-1).*r12.*(ud2+(-1).*ud6+x02+(-1).*x06)).*(abs(vw1).^2+ ...
  abs(vw2).^2).*(abs(ud1+(-1).*ud5+x01+(-1).*x05).^2+abs(ud2+(-1).* ...
  ud6+x02+(-1).*x06).^2).^(-1)+(-32).*A.*rho.*vw1.*(4.*vw1.^2+4.* ...
  vw2.^2).^(-3).*(vw1.*(r11.*(ud1+(-1).*ud5+x01+(-1).*x05)+r12.*( ...
  ud2+(-1).*ud6+x02+(-1).*x06))+vw2.*(r21.*(ud1+(-1).*ud5+x01+(-1).* ...
  x05)+r22.*(ud2+(-1).*ud6+x02+(-1).*x06))).^2.*abs(vw1).*(abs(ud1+( ...
  -1).*ud5+x01+(-1).*x05).^2+abs(ud2+(-1).*ud6+x02+(-1).*x06).^2).^( ...
  -1).*sign(vw1),384.*A.*rho.*vw1.*vw2.*(4.*vw1.^2+4.*vw2.^2).^(-4) ...
  .*(vw1.*(r11.*(ud1+(-1).*ud5+x01+(-1).*x05)+r12.*(ud2+(-1).*ud6+ ...
  x02+(-1).*x06))+vw2.*(r21.*(ud1+(-1).*ud5+x01+(-1).*x05)+r22.*( ...
  ud2+(-1).*ud6+x02+(-1).*x06))).^2.*(abs(vw1).^2+abs(vw2).^2).*( ...
  abs(ud1+(-1).*ud5+x01+(-1).*x05).^2+abs(ud2+(-1).*ud6+x02+(-1).* ...
  x06).^2).^(-1)+32.*A.*rho.*vw1.*(4.*vw1.^2+4.*vw2.^2).^(-3).*( ...
  vw1.*(r11.*(ud1+(-1).*ud5+x01+(-1).*x05)+r12.*(ud2+(-1).*ud6+x02+( ...
  -1).*x06))+vw2.*(r21.*(ud1+(-1).*ud5+x01+(-1).*x05)+r22.*(ud2+(-1) ...
  .*ud6+x02+(-1).*x06))).*((-1).*r21.*(ud1+(-1).*ud5+x01+(-1).*x05)+ ...
  (-1).*r22.*(ud2+(-1).*ud6+x02+(-1).*x06)).*(abs(vw1).^2+abs(vw2) ...
  .^2).*(abs(ud1+(-1).*ud5+x01+(-1).*x05).^2+abs(ud2+(-1).*ud6+x02+( ...
  -1).*x06).^2).^(-1)+(-32).*A.*rho.*vw1.*(4.*vw1.^2+4.*vw2.^2).^( ...
  -3).*(vw1.*(r11.*(ud1+(-1).*ud5+x01+(-1).*x05)+r12.*(ud2+(-1).* ...
  ud6+x02+(-1).*x06))+vw2.*(r21.*(ud1+(-1).*ud5+x01+(-1).*x05)+r22.* ...
  (ud2+(-1).*ud6+x02+(-1).*x06))).^2.*abs(vw2).*(abs(ud1+(-1).*ud5+ ...
  x01+(-1).*x05).^2+abs(ud2+(-1).*ud6+x02+(-1).*x06).^2).^(-1).* ...
  sign(vw2);384.*A.*rho.*vw1.*vw2.*(4.*vw1.^2+4.*vw2.^2).^(-4).*( ...
  vw1.*(r11.*(ud1+(-1).*ud5+x01+(-1).*x05)+r12.*(ud2+(-1).*ud6+x02+( ...
  -1).*x06))+vw2.*(r21.*(ud1+(-1).*ud5+x01+(-1).*x05)+r22.*(ud2+(-1) ...
  .*ud6+x02+(-1).*x06))).^2.*(abs(vw1).^2+abs(vw2).^2).*(abs(ud1+( ...
  -1).*ud5+x01+(-1).*x05).^2+abs(ud2+(-1).*ud6+x02+(-1).*x06).^2).^( ...
  -1)+32.*A.*rho.*vw2.*(4.*vw1.^2+4.*vw2.^2).^(-3).*(vw1.*(r11.*( ...
  ud1+(-1).*ud5+x01+(-1).*x05)+r12.*(ud2+(-1).*ud6+x02+(-1).*x06))+ ...
  vw2.*(r21.*(ud1+(-1).*ud5+x01+(-1).*x05)+r22.*(ud2+(-1).*ud6+x02+( ...
  -1).*x06))).*((-1).*r11.*(ud1+(-1).*ud5+x01+(-1).*x05)+(-1).*r12.* ...
  (ud2+(-1).*ud6+x02+(-1).*x06)).*(abs(vw1).^2+abs(vw2).^2).*(abs( ...
  ud1+(-1).*ud5+x01+(-1).*x05).^2+abs(ud2+(-1).*ud6+x02+(-1).*x06) ...
  .^2).^(-1)+(-32).*A.*rho.*vw2.*(4.*vw1.^2+4.*vw2.^2).^(-3).*(vw1.* ...
  (r11.*(ud1+(-1).*ud5+x01+(-1).*x05)+r12.*(ud2+(-1).*ud6+x02+(-1).* ...
  x06))+vw2.*(r21.*(ud1+(-1).*ud5+x01+(-1).*x05)+r22.*(ud2+(-1).* ...
  ud6+x02+(-1).*x06))).^2.*abs(vw1).*(abs(ud1+(-1).*ud5+x01+(-1).* ...
  x05).^2+abs(ud2+(-1).*ud6+x02+(-1).*x06).^2).^(-1).*sign(vw1), ...
  384.*A.*rho.*vw2.^2.*(4.*vw1.^2+4.*vw2.^2).^(-4).*(vw1.*(r11.*( ...
  ud1+(-1).*ud5+x01+(-1).*x05)+r12.*(ud2+(-1).*ud6+x02+(-1).*x06))+ ...
  vw2.*(r21.*(ud1+(-1).*ud5+x01+(-1).*x05)+r22.*(ud2+(-1).*ud6+x02+( ...
  -1).*x06))).^2.*(abs(vw1).^2+abs(vw2).^2).*(abs(ud1+(-1).*ud5+x01+ ...
  (-1).*x05).^2+abs(ud2+(-1).*ud6+x02+(-1).*x06).^2).^(-1)+(-16).* ...
  A.*rho.*(4.*vw1.^2+4.*vw2.^2).^(-3).*(vw1.*(r11.*(ud1+(-1).*ud5+ ...
  x01+(-1).*x05)+r12.*(ud2+(-1).*ud6+x02+(-1).*x06))+vw2.*(r21.*( ...
  ud1+(-1).*ud5+x01+(-1).*x05)+r22.*(ud2+(-1).*ud6+x02+(-1).*x06))) ...
  .^2.*(abs(vw1).^2+abs(vw2).^2).*(abs(ud1+(-1).*ud5+x01+(-1).*x05) ...
  .^2+abs(ud2+(-1).*ud6+x02+(-1).*x06).^2).^(-1)+32.*A.*rho.*vw2.*( ...
  4.*vw1.^2+4.*vw2.^2).^(-3).*(vw1.*(r11.*(ud1+(-1).*ud5+x01+(-1).* ...
  x05)+r12.*(ud2+(-1).*ud6+x02+(-1).*x06))+vw2.*(r21.*(ud1+(-1).* ...
  ud5+x01+(-1).*x05)+r22.*(ud2+(-1).*ud6+x02+(-1).*x06))).*((-1).* ...
  r21.*(ud1+(-1).*ud5+x01+(-1).*x05)+(-1).*r22.*(ud2+(-1).*ud6+x02+( ...
  -1).*x06)).*(abs(vw1).^2+abs(vw2).^2).*(abs(ud1+(-1).*ud5+x01+(-1) ...
  .*x05).^2+abs(ud2+(-1).*ud6+x02+(-1).*x06).^2).^(-1)+(-32).*A.* ...
  rho.*vw2.*(4.*vw1.^2+4.*vw2.^2).^(-3).*(vw1.*(r11.*(ud1+(-1).*ud5+ ...
  x01+(-1).*x05)+r12.*(ud2+(-1).*ud6+x02+(-1).*x06))+vw2.*(r21.*( ...
  ud1+(-1).*ud5+x01+(-1).*x05)+r22.*(ud2+(-1).*ud6+x02+(-1).*x06))) ...
  .^2.*abs(vw2).*(abs(ud1+(-1).*ud5+x01+(-1).*x05).^2+abs(ud2+(-1).* ...
  ud6+x02+(-1).*x06).^2).^(-1).*sign(vw2),0,0,384.*A.*rho.*vw1.* ...
  vw2.*(4.*vw1.^2+4.*vw2.^2).^(-4).*(vw1.*(r11.*(ud1+(-1).*ud5+x01+( ...
  -1).*x05)+r12.*(ud2+(-1).*ud6+x02+(-1).*x06))+vw2.*(r21.*(ud1+(-1) ...
  .*ud5+x01+(-1).*x05)+r22.*(ud2+(-1).*ud6+x02+(-1).*x06))).^2.*( ...
  abs(vw1).^2+abs(vw2).^2).*(abs(ud1+(-1).*ud5+x01+(-1).*x05).^2+ ...
  abs(ud2+(-1).*ud6+x02+(-1).*x06).^2).^(-1)+32.*A.*rho.*vw2.*(4.* ...
  vw1.^2+4.*vw2.^2).^(-3).*(vw1.*(r11.*(ud1+(-1).*ud5+x01+(-1).*x05) ...
  +r12.*(ud2+(-1).*ud6+x02+(-1).*x06))+vw2.*(r21.*(ud1+(-1).*ud5+ ...
  x01+(-1).*x05)+r22.*(ud2+(-1).*ud6+x02+(-1).*x06))).*((-1).*r11.*( ...
  ud1+(-1).*ud5+x01+(-1).*x05)+(-1).*r12.*(ud2+(-1).*ud6+x02+(-1).* ...
  x06)).*(abs(vw1).^2+abs(vw2).^2).*(abs(ud1+(-1).*ud5+x01+(-1).* ...
  x05).^2+abs(ud2+(-1).*ud6+x02+(-1).*x06).^2).^(-1)+(-32).*A.*rho.* ...
  vw2.*(4.*vw1.^2+4.*vw2.^2).^(-3).*(vw1.*(r11.*(ud1+(-1).*ud5+x01+( ...
  -1).*x05)+r12.*(ud2+(-1).*ud6+x02+(-1).*x06))+vw2.*(r21.*(ud1+(-1) ...
  .*ud5+x01+(-1).*x05)+r22.*(ud2+(-1).*ud6+x02+(-1).*x06))).^2.*abs( ...
  vw1).*(abs(ud1+(-1).*ud5+x01+(-1).*x05).^2+abs(ud2+(-1).*ud6+x02+( ...
  -1).*x06).^2).^(-1).*sign(vw1),384.*A.*rho.*vw2.^2.*(4.*vw1.^2+4.* ...
  vw2.^2).^(-4).*(vw1.*(r11.*(ud1+(-1).*ud5+x01+(-1).*x05)+r12.*( ...
  ud2+(-1).*ud6+x02+(-1).*x06))+vw2.*(r21.*(ud1+(-1).*ud5+x01+(-1).* ...
  x05)+r22.*(ud2+(-1).*ud6+x02+(-1).*x06))).^2.*(abs(vw1).^2+abs( ...
  vw2).^2).*(abs(ud1+(-1).*ud5+x01+(-1).*x05).^2+abs(ud2+(-1).*ud6+ ...
  x02+(-1).*x06).^2).^(-1)+(-16).*A.*rho.*(4.*vw1.^2+4.*vw2.^2).^( ...
  -3).*(vw1.*(r11.*(ud1+(-1).*ud5+x01+(-1).*x05)+r12.*(ud2+(-1).* ...
  ud6+x02+(-1).*x06))+vw2.*(r21.*(ud1+(-1).*ud5+x01+(-1).*x05)+r22.* ...
  (ud2+(-1).*ud6+x02+(-1).*x06))).^2.*(abs(vw1).^2+abs(vw2).^2).*( ...
  abs(ud1+(-1).*ud5+x01+(-1).*x05).^2+abs(ud2+(-1).*ud6+x02+(-1).* ...
  x06).^2).^(-1)+32.*A.*rho.*vw2.*(4.*vw1.^2+4.*vw2.^2).^(-3).*( ...
  vw1.*(r11.*(ud1+(-1).*ud5+x01+(-1).*x05)+r12.*(ud2+(-1).*ud6+x02+( ...
  -1).*x06))+vw2.*(r21.*(ud1+(-1).*ud5+x01+(-1).*x05)+r22.*(ud2+(-1) ...
  .*ud6+x02+(-1).*x06))).*((-1).*r21.*(ud1+(-1).*ud5+x01+(-1).*x05)+ ...
  (-1).*r22.*(ud2+(-1).*ud6+x02+(-1).*x06)).*(abs(vw1).^2+abs(vw2) ...
  .^2).*(abs(ud1+(-1).*ud5+x01+(-1).*x05).^2+abs(ud2+(-1).*ud6+x02+( ...
  -1).*x06).^2).^(-1)+(-32).*A.*rho.*vw2.*(4.*vw1.^2+4.*vw2.^2).^( ...
  -3).*(vw1.*(r11.*(ud1+(-1).*ud5+x01+(-1).*x05)+r12.*(ud2+(-1).* ...
  ud6+x02+(-1).*x06))+vw2.*(r21.*(ud1+(-1).*ud5+x01+(-1).*x05)+r22.* ...
  (ud2+(-1).*ud6+x02+(-1).*x06))).^2.*abs(vw2).*(abs(ud1+(-1).*ud5+ ...
  x01+(-1).*x05).^2+abs(ud2+(-1).*ud6+x02+(-1).*x06).^2).^(-1).* ...
  sign(vw2);0,0,0,0,0,0;0,0,0,0,0,0;384.*A.*rho.*vw1.^2.*(4.*vw1.^2+ ...
  4.*vw2.^2).^(-4).*(vw1.*(r11.*(ud1+(-1).*ud5+x01+(-1).*x05)+r12.*( ...
  ud2+(-1).*ud6+x02+(-1).*x06))+vw2.*(r21.*(ud1+(-1).*ud5+x01+(-1).* ...
  x05)+r22.*(ud2+(-1).*ud6+x02+(-1).*x06))).^2.*(abs(vw1).^2+abs( ...
  vw2).^2).*(abs(ud1+(-1).*ud5+x01+(-1).*x05).^2+abs(ud2+(-1).*ud6+ ...
  x02+(-1).*x06).^2).^(-1)+(-16).*A.*rho.*(4.*vw1.^2+4.*vw2.^2).^( ...
  -3).*(vw1.*(r11.*(ud1+(-1).*ud5+x01+(-1).*x05)+r12.*(ud2+(-1).* ...
  ud6+x02+(-1).*x06))+vw2.*(r21.*(ud1+(-1).*ud5+x01+(-1).*x05)+r22.* ...
  (ud2+(-1).*ud6+x02+(-1).*x06))).^2.*(abs(vw1).^2+abs(vw2).^2).*( ...
  abs(ud1+(-1).*ud5+x01+(-1).*x05).^2+abs(ud2+(-1).*ud6+x02+(-1).* ...
  x06).^2).^(-1)+32.*A.*rho.*vw1.*(4.*vw1.^2+4.*vw2.^2).^(-3).*( ...
  vw1.*(r11.*(ud1+(-1).*ud5+x01+(-1).*x05)+r12.*(ud2+(-1).*ud6+x02+( ...
  -1).*x06))+vw2.*(r21.*(ud1+(-1).*ud5+x01+(-1).*x05)+r22.*(ud2+(-1) ...
  .*ud6+x02+(-1).*x06))).*((-1).*r11.*(ud1+(-1).*ud5+x01+(-1).*x05)+ ...
  (-1).*r12.*(ud2+(-1).*ud6+x02+(-1).*x06)).*(abs(vw1).^2+abs(vw2) ...
  .^2).*(abs(ud1+(-1).*ud5+x01+(-1).*x05).^2+abs(ud2+(-1).*ud6+x02+( ...
  -1).*x06).^2).^(-1)+(-32).*A.*rho.*vw1.*(4.*vw1.^2+4.*vw2.^2).^( ...
  -3).*(vw1.*(r11.*(ud1+(-1).*ud5+x01+(-1).*x05)+r12.*(ud2+(-1).* ...
  ud6+x02+(-1).*x06))+vw2.*(r21.*(ud1+(-1).*ud5+x01+(-1).*x05)+r22.* ...
  (ud2+(-1).*ud6+x02+(-1).*x06))).^2.*abs(vw1).*(abs(ud1+(-1).*ud5+ ...
  x01+(-1).*x05).^2+abs(ud2+(-1).*ud6+x02+(-1).*x06).^2).^(-1).* ...
  sign(vw1),384.*A.*rho.*vw1.*vw2.*(4.*vw1.^2+4.*vw2.^2).^(-4).*( ...
  vw1.*(r11.*(ud1+(-1).*ud5+x01+(-1).*x05)+r12.*(ud2+(-1).*ud6+x02+( ...
  -1).*x06))+vw2.*(r21.*(ud1+(-1).*ud5+x01+(-1).*x05)+r22.*(ud2+(-1) ...
  .*ud6+x02+(-1).*x06))).^2.*(abs(vw1).^2+abs(vw2).^2).*(abs(ud1+( ...
  -1).*ud5+x01+(-1).*x05).^2+abs(ud2+(-1).*ud6+x02+(-1).*x06).^2).^( ...
  -1)+32.*A.*rho.*vw1.*(4.*vw1.^2+4.*vw2.^2).^(-3).*(vw1.*(r11.*( ...
  ud1+(-1).*ud5+x01+(-1).*x05)+r12.*(ud2+(-1).*ud6+x02+(-1).*x06))+ ...
  vw2.*(r21.*(ud1+(-1).*ud5+x01+(-1).*x05)+r22.*(ud2+(-1).*ud6+x02+( ...
  -1).*x06))).*((-1).*r21.*(ud1+(-1).*ud5+x01+(-1).*x05)+(-1).*r22.* ...
  (ud2+(-1).*ud6+x02+(-1).*x06)).*(abs(vw1).^2+abs(vw2).^2).*(abs( ...
  ud1+(-1).*ud5+x01+(-1).*x05).^2+abs(ud2+(-1).*ud6+x02+(-1).*x06) ...
  .^2).^(-1)+(-32).*A.*rho.*vw1.*(4.*vw1.^2+4.*vw2.^2).^(-3).*(vw1.* ...
  (r11.*(ud1+(-1).*ud5+x01+(-1).*x05)+r12.*(ud2+(-1).*ud6+x02+(-1).* ...
  x06))+vw2.*(r21.*(ud1+(-1).*ud5+x01+(-1).*x05)+r22.*(ud2+(-1).* ...
  ud6+x02+(-1).*x06))).^2.*abs(vw2).*(abs(ud1+(-1).*ud5+x01+(-1).* ...
  x05).^2+abs(ud2+(-1).*ud6+x02+(-1).*x06).^2).^(-1).*sign(vw2),0,0, ...
  384.*A.*rho.*vw1.^2.*(4.*vw1.^2+4.*vw2.^2).^(-4).*(vw1.*(r11.*( ...
  ud1+(-1).*ud5+x01+(-1).*x05)+r12.*(ud2+(-1).*ud6+x02+(-1).*x06))+ ...
  vw2.*(r21.*(ud1+(-1).*ud5+x01+(-1).*x05)+r22.*(ud2+(-1).*ud6+x02+( ...
  -1).*x06))).^2.*(abs(vw1).^2+abs(vw2).^2).*(abs(ud1+(-1).*ud5+x01+ ...
  (-1).*x05).^2+abs(ud2+(-1).*ud6+x02+(-1).*x06).^2).^(-1)+(-16).* ...
  A.*rho.*(4.*vw1.^2+4.*vw2.^2).^(-3).*(vw1.*(r11.*(ud1+(-1).*ud5+ ...
  x01+(-1).*x05)+r12.*(ud2+(-1).*ud6+x02+(-1).*x06))+vw2.*(r21.*( ...
  ud1+(-1).*ud5+x01+(-1).*x05)+r22.*(ud2+(-1).*ud6+x02+(-1).*x06))) ...
  .^2.*(abs(vw1).^2+abs(vw2).^2).*(abs(ud1+(-1).*ud5+x01+(-1).*x05) ...
  .^2+abs(ud2+(-1).*ud6+x02+(-1).*x06).^2).^(-1)+32.*A.*rho.*vw1.*( ...
  4.*vw1.^2+4.*vw2.^2).^(-3).*(vw1.*(r11.*(ud1+(-1).*ud5+x01+(-1).* ...
  x05)+r12.*(ud2+(-1).*ud6+x02+(-1).*x06))+vw2.*(r21.*(ud1+(-1).* ...
  ud5+x01+(-1).*x05)+r22.*(ud2+(-1).*ud6+x02+(-1).*x06))).*((-1).* ...
  r11.*(ud1+(-1).*ud5+x01+(-1).*x05)+(-1).*r12.*(ud2+(-1).*ud6+x02+( ...
  -1).*x06)).*(abs(vw1).^2+abs(vw2).^2).*(abs(ud1+(-1).*ud5+x01+(-1) ...
  .*x05).^2+abs(ud2+(-1).*ud6+x02+(-1).*x06).^2).^(-1)+(-32).*A.* ...
  rho.*vw1.*(4.*vw1.^2+4.*vw2.^2).^(-3).*(vw1.*(r11.*(ud1+(-1).*ud5+ ...
  x01+(-1).*x05)+r12.*(ud2+(-1).*ud6+x02+(-1).*x06))+vw2.*(r21.*( ...
  ud1+(-1).*ud5+x01+(-1).*x05)+r22.*(ud2+(-1).*ud6+x02+(-1).*x06))) ...
  .^2.*abs(vw1).*(abs(ud1+(-1).*ud5+x01+(-1).*x05).^2+abs(ud2+(-1).* ...
  ud6+x02+(-1).*x06).^2).^(-1).*sign(vw1),384.*A.*rho.*vw1.*vw2.*( ...
  4.*vw1.^2+4.*vw2.^2).^(-4).*(vw1.*(r11.*(ud1+(-1).*ud5+x01+(-1).* ...
  x05)+r12.*(ud2+(-1).*ud6+x02+(-1).*x06))+vw2.*(r21.*(ud1+(-1).* ...
  ud5+x01+(-1).*x05)+r22.*(ud2+(-1).*ud6+x02+(-1).*x06))).^2.*(abs( ...
  vw1).^2+abs(vw2).^2).*(abs(ud1+(-1).*ud5+x01+(-1).*x05).^2+abs( ...
  ud2+(-1).*ud6+x02+(-1).*x06).^2).^(-1)+32.*A.*rho.*vw1.*(4.* ...
  vw1.^2+4.*vw2.^2).^(-3).*(vw1.*(r11.*(ud1+(-1).*ud5+x01+(-1).*x05) ...
  +r12.*(ud2+(-1).*ud6+x02+(-1).*x06))+vw2.*(r21.*(ud1+(-1).*ud5+ ...
  x01+(-1).*x05)+r22.*(ud2+(-1).*ud6+x02+(-1).*x06))).*((-1).*r21.*( ...
  ud1+(-1).*ud5+x01+(-1).*x05)+(-1).*r22.*(ud2+(-1).*ud6+x02+(-1).* ...
  x06)).*(abs(vw1).^2+abs(vw2).^2).*(abs(ud1+(-1).*ud5+x01+(-1).* ...
  x05).^2+abs(ud2+(-1).*ud6+x02+(-1).*x06).^2).^(-1)+(-32).*A.*rho.* ...
  vw1.*(4.*vw1.^2+4.*vw2.^2).^(-3).*(vw1.*(r11.*(ud1+(-1).*ud5+x01+( ...
  -1).*x05)+r12.*(ud2+(-1).*ud6+x02+(-1).*x06))+vw2.*(r21.*(ud1+(-1) ...
  .*ud5+x01+(-1).*x05)+r22.*(ud2+(-1).*ud6+x02+(-1).*x06))).^2.*abs( ...
  vw2).*(abs(ud1+(-1).*ud5+x01+(-1).*x05).^2+abs(ud2+(-1).*ud6+x02+( ...
  -1).*x06).^2).^(-1).*sign(vw2);384.*A.*rho.*vw1.*vw2.*(4.*vw1.^2+ ...
  4.*vw2.^2).^(-4).*(vw1.*(r11.*(ud1+(-1).*ud5+x01+(-1).*x05)+r12.*( ...
  ud2+(-1).*ud6+x02+(-1).*x06))+vw2.*(r21.*(ud1+(-1).*ud5+x01+(-1).* ...
  x05)+r22.*(ud2+(-1).*ud6+x02+(-1).*x06))).^2.*(abs(vw1).^2+abs( ...
  vw2).^2).*(abs(ud1+(-1).*ud5+x01+(-1).*x05).^2+abs(ud2+(-1).*ud6+ ...
  x02+(-1).*x06).^2).^(-1)+32.*A.*rho.*vw2.*(4.*vw1.^2+4.*vw2.^2).^( ...
  -3).*(vw1.*(r11.*(ud1+(-1).*ud5+x01+(-1).*x05)+r12.*(ud2+(-1).* ...
  ud6+x02+(-1).*x06))+vw2.*(r21.*(ud1+(-1).*ud5+x01+(-1).*x05)+r22.* ...
  (ud2+(-1).*ud6+x02+(-1).*x06))).*((-1).*r11.*(ud1+(-1).*ud5+x01+( ...
  -1).*x05)+(-1).*r12.*(ud2+(-1).*ud6+x02+(-1).*x06)).*(abs(vw1).^2+ ...
  abs(vw2).^2).*(abs(ud1+(-1).*ud5+x01+(-1).*x05).^2+abs(ud2+(-1).* ...
  ud6+x02+(-1).*x06).^2).^(-1)+(-32).*A.*rho.*vw2.*(4.*vw1.^2+4.* ...
  vw2.^2).^(-3).*(vw1.*(r11.*(ud1+(-1).*ud5+x01+(-1).*x05)+r12.*( ...
  ud2+(-1).*ud6+x02+(-1).*x06))+vw2.*(r21.*(ud1+(-1).*ud5+x01+(-1).* ...
  x05)+r22.*(ud2+(-1).*ud6+x02+(-1).*x06))).^2.*abs(vw1).*(abs(ud1+( ...
  -1).*ud5+x01+(-1).*x05).^2+abs(ud2+(-1).*ud6+x02+(-1).*x06).^2).^( ...
  -1).*sign(vw1),384.*A.*rho.*vw2.^2.*(4.*vw1.^2+4.*vw2.^2).^(-4).*( ...
  vw1.*(r11.*(ud1+(-1).*ud5+x01+(-1).*x05)+r12.*(ud2+(-1).*ud6+x02+( ...
  -1).*x06))+vw2.*(r21.*(ud1+(-1).*ud5+x01+(-1).*x05)+r22.*(ud2+(-1) ...
  .*ud6+x02+(-1).*x06))).^2.*(abs(vw1).^2+abs(vw2).^2).*(abs(ud1+( ...
  -1).*ud5+x01+(-1).*x05).^2+abs(ud2+(-1).*ud6+x02+(-1).*x06).^2).^( ...
  -1)+(-16).*A.*rho.*(4.*vw1.^2+4.*vw2.^2).^(-3).*(vw1.*(r11.*(ud1+( ...
  -1).*ud5+x01+(-1).*x05)+r12.*(ud2+(-1).*ud6+x02+(-1).*x06))+vw2.*( ...
  r21.*(ud1+(-1).*ud5+x01+(-1).*x05)+r22.*(ud2+(-1).*ud6+x02+(-1).* ...
  x06))).^2.*(abs(vw1).^2+abs(vw2).^2).*(abs(ud1+(-1).*ud5+x01+(-1) ...
  .*x05).^2+abs(ud2+(-1).*ud6+x02+(-1).*x06).^2).^(-1)+32.*A.*rho.* ...
  vw2.*(4.*vw1.^2+4.*vw2.^2).^(-3).*(vw1.*(r11.*(ud1+(-1).*ud5+x01+( ...
  -1).*x05)+r12.*(ud2+(-1).*ud6+x02+(-1).*x06))+vw2.*(r21.*(ud1+(-1) ...
  .*ud5+x01+(-1).*x05)+r22.*(ud2+(-1).*ud6+x02+(-1).*x06))).*((-1).* ...
  r21.*(ud1+(-1).*ud5+x01+(-1).*x05)+(-1).*r22.*(ud2+(-1).*ud6+x02+( ...
  -1).*x06)).*(abs(vw1).^2+abs(vw2).^2).*(abs(ud1+(-1).*ud5+x01+(-1) ...
  .*x05).^2+abs(ud2+(-1).*ud6+x02+(-1).*x06).^2).^(-1)+(-32).*A.* ...
  rho.*vw2.*(4.*vw1.^2+4.*vw2.^2).^(-3).*(vw1.*(r11.*(ud1+(-1).*ud5+ ...
  x01+(-1).*x05)+r12.*(ud2+(-1).*ud6+x02+(-1).*x06))+vw2.*(r21.*( ...
  ud1+(-1).*ud5+x01+(-1).*x05)+r22.*(ud2+(-1).*ud6+x02+(-1).*x06))) ...
  .^2.*abs(vw2).*(abs(ud1+(-1).*ud5+x01+(-1).*x05).^2+abs(ud2+(-1).* ...
  ud6+x02+(-1).*x06).^2).^(-1).*sign(vw2),0,0,384.*A.*rho.*vw1.* ...
  vw2.*(4.*vw1.^2+4.*vw2.^2).^(-4).*(vw1.*(r11.*(ud1+(-1).*ud5+x01+( ...
  -1).*x05)+r12.*(ud2+(-1).*ud6+x02+(-1).*x06))+vw2.*(r21.*(ud1+(-1) ...
  .*ud5+x01+(-1).*x05)+r22.*(ud2+(-1).*ud6+x02+(-1).*x06))).^2.*( ...
  abs(vw1).^2+abs(vw2).^2).*(abs(ud1+(-1).*ud5+x01+(-1).*x05).^2+ ...
  abs(ud2+(-1).*ud6+x02+(-1).*x06).^2).^(-1)+32.*A.*rho.*vw2.*(4.* ...
  vw1.^2+4.*vw2.^2).^(-3).*(vw1.*(r11.*(ud1+(-1).*ud5+x01+(-1).*x05) ...
  +r12.*(ud2+(-1).*ud6+x02+(-1).*x06))+vw2.*(r21.*(ud1+(-1).*ud5+ ...
  x01+(-1).*x05)+r22.*(ud2+(-1).*ud6+x02+(-1).*x06))).*((-1).*r11.*( ...
  ud1+(-1).*ud5+x01+(-1).*x05)+(-1).*r12.*(ud2+(-1).*ud6+x02+(-1).* ...
  x06)).*(abs(vw1).^2+abs(vw2).^2).*(abs(ud1+(-1).*ud5+x01+(-1).* ...
  x05).^2+abs(ud2+(-1).*ud6+x02+(-1).*x06).^2).^(-1)+(-32).*A.*rho.* ...
  vw2.*(4.*vw1.^2+4.*vw2.^2).^(-3).*(vw1.*(r11.*(ud1+(-1).*ud5+x01+( ...
  -1).*x05)+r12.*(ud2+(-1).*ud6+x02+(-1).*x06))+vw2.*(r21.*(ud1+(-1) ...
  .*ud5+x01+(-1).*x05)+r22.*(ud2+(-1).*ud6+x02+(-1).*x06))).^2.*abs( ...
  vw1).*(abs(ud1+(-1).*ud5+x01+(-1).*x05).^2+abs(ud2+(-1).*ud6+x02+( ...
  -1).*x06).^2).^(-1).*sign(vw1),384.*A.*rho.*vw2.^2.*(4.*vw1.^2+4.* ...
  vw2.^2).^(-4).*(vw1.*(r11.*(ud1+(-1).*ud5+x01+(-1).*x05)+r12.*( ...
  ud2+(-1).*ud6+x02+(-1).*x06))+vw2.*(r21.*(ud1+(-1).*ud5+x01+(-1).* ...
  x05)+r22.*(ud2+(-1).*ud6+x02+(-1).*x06))).^2.*(abs(vw1).^2+abs( ...
  vw2).^2).*(abs(ud1+(-1).*ud5+x01+(-1).*x05).^2+abs(ud2+(-1).*ud6+ ...
  x02+(-1).*x06).^2).^(-1)+(-16).*A.*rho.*(4.*vw1.^2+4.*vw2.^2).^( ...
  -3).*(vw1.*(r11.*(ud1+(-1).*ud5+x01+(-1).*x05)+r12.*(ud2+(-1).* ...
  ud6+x02+(-1).*x06))+vw2.*(r21.*(ud1+(-1).*ud5+x01+(-1).*x05)+r22.* ...
  (ud2+(-1).*ud6+x02+(-1).*x06))).^2.*(abs(vw1).^2+abs(vw2).^2).*( ...
  abs(ud1+(-1).*ud5+x01+(-1).*x05).^2+abs(ud2+(-1).*ud6+x02+(-1).* ...
  x06).^2).^(-1)+32.*A.*rho.*vw2.*(4.*vw1.^2+4.*vw2.^2).^(-3).*( ...
  vw1.*(r11.*(ud1+(-1).*ud5+x01+(-1).*x05)+r12.*(ud2+(-1).*ud6+x02+( ...
  -1).*x06))+vw2.*(r21.*(ud1+(-1).*ud5+x01+(-1).*x05)+r22.*(ud2+(-1) ...
  .*ud6+x02+(-1).*x06))).*((-1).*r21.*(ud1+(-1).*ud5+x01+(-1).*x05)+ ...
  (-1).*r22.*(ud2+(-1).*ud6+x02+(-1).*x06)).*(abs(vw1).^2+abs(vw2) ...
  .^2).*(abs(ud1+(-1).*ud5+x01+(-1).*x05).^2+abs(ud2+(-1).*ud6+x02+( ...
  -1).*x06).^2).^(-1)+(-32).*A.*rho.*vw2.*(4.*vw1.^2+4.*vw2.^2).^( ...
  -3).*(vw1.*(r11.*(ud1+(-1).*ud5+x01+(-1).*x05)+r12.*(ud2+(-1).* ...
  ud6+x02+(-1).*x06))+vw2.*(r21.*(ud1+(-1).*ud5+x01+(-1).*x05)+r22.* ...
  (ud2+(-1).*ud6+x02+(-1).*x06))).^2.*abs(vw2).*(abs(ud1+(-1).*ud5+ ...
  x01+(-1).*x05).^2+abs(ud2+(-1).*ud6+x02+(-1).*x06).^2).^(-1).* ...
  sign(vw2)];


end
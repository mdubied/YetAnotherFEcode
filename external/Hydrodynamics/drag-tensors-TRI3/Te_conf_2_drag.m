function T = Te_conf_2_drag(A,rho,vw1,vw2,r11,r12,r21,r22,x01,x02,x03,x04,x05,x06,ud1,ud2,ud3,ud4,ud5,ud6)
    T = [0,0,32.*A.*rho.*vw1.*(4.*vw1.^2+4.*vw2.^2).^(-3).*(vw1.*(r11.*(( ...
  -1).*ud3+ud5+(-1).*x03+x05)+r12.*((-1).*ud4+ud6+(-1).*x04+x06))+ ...
  vw2.*(r21.*((-1).*ud3+ud5+(-1).*x03+x05)+r22.*((-1).*ud4+ud6+(-1) ...
  .*x04+x06))).^2.*(abs(vw1).^2+abs(vw2).^2).*(abs(ud3+(-1).*ud5+ ...
  x03+(-1).*x05).^2+abs(ud4+(-1).*ud6+x04+(-1).*x06).^2).^(-1),32.* ...
  A.*rho.*vw2.*(4.*vw1.^2+4.*vw2.^2).^(-3).*(vw1.*(r11.*((-1).*ud3+ ...
  ud5+(-1).*x03+x05)+r12.*((-1).*ud4+ud6+(-1).*x04+x06))+vw2.*(r21.* ...
  ((-1).*ud3+ud5+(-1).*x03+x05)+r22.*((-1).*ud4+ud6+(-1).*x04+x06))) ...
  .^2.*(abs(vw1).^2+abs(vw2).^2).*(abs(ud3+(-1).*ud5+x03+(-1).*x05) ...
  .^2+abs(ud4+(-1).*ud6+x04+(-1).*x06).^2).^(-1),32.*A.*rho.*vw1.*( ...
  4.*vw1.^2+4.*vw2.^2).^(-3).*(vw1.*(r11.*((-1).*ud3+ud5+(-1).*x03+ ...
  x05)+r12.*((-1).*ud4+ud6+(-1).*x04+x06))+vw2.*(r21.*((-1).*ud3+ ...
  ud5+(-1).*x03+x05)+r22.*((-1).*ud4+ud6+(-1).*x04+x06))).^2.*(abs( ...
  vw1).^2+abs(vw2).^2).*(abs(ud3+(-1).*ud5+x03+(-1).*x05).^2+abs( ...
  ud4+(-1).*ud6+x04+(-1).*x06).^2).^(-1),32.*A.*rho.*vw2.*(4.* ...
  vw1.^2+4.*vw2.^2).^(-3).*(vw1.*(r11.*((-1).*ud3+ud5+(-1).*x03+x05) ...
  +r12.*((-1).*ud4+ud6+(-1).*x04+x06))+vw2.*(r21.*((-1).*ud3+ud5+( ...
  -1).*x03+x05)+r22.*((-1).*ud4+ud6+(-1).*x04+x06))).^2.*(abs(vw1) ...
  .^2+abs(vw2).^2).*(abs(ud3+(-1).*ud5+x03+(-1).*x05).^2+abs(ud4+( ...
  -1).*ud6+x04+(-1).*x06).^2).^(-1)];
end 
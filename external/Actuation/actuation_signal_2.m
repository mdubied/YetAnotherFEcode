% actuation_signal_2
%
% Syntax: actuation_signal_2(k,t,p)
%
% Description: function producing the actuation signal 2 described in the
% paper. The actuation signal is defined as df_actu/dq
%
% Last modified: 12/10/2023, Mathieu Dubied, ETH Zürich
function actu = actuation_signal_2(k,t,p)
    p1 = p(1);
    p2 = p(2);
    p3 = p(3);
    actu =  p1*k/2*(-0.2*sin(t*2*pi)) + ...
            p2*k/2*(-0.2*sin(t*2*pi*3)) + ...
            p3*k/2*(-0.2*sin(t*2*pi*5));
end
% actuation_signal_1
%
% Syntax: actuation_signal_1(k,t,p1)
%
% Description: function producing the actuation signal 1 described in the
% paper
%
% Last modified: 10/10/2023, Mathieu Dubied, ETH Zürich
function actu = actuation_signal_1(k,t,p1)
    actu = p1*k/2*(-0.2*sin(t*2*pi));
end
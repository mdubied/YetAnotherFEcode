% DpROM, reference: 
% Marconi et al. (2021). "A higher-order parametric nonlinear 
% reduced-order model for imperfect structures using Neumann 
% expansion". Nonlinear Dynamics. 
% https://doi.org/10.1007/s11071-021-06496-y

function A2 = A2_fun(self, th)
    if self.nDim == 2
        A2 = -[th(1) th(3) 0     0;
              0     0     th(2) th(4);
              th(2) th(4) th(1) th(3)];
    elseif self.nDim == 3
        A2 = -[ ...
        th(1),  th(4),  th(7),      0,      0,      0,      0,      0,      0;
            0,      0,      0,  th(2),  th(5),  th(8),      0,      0,      0;
            0,      0,      0,      0,      0,      0,  th(3),  th(6),  th(9);
        th(2),  th(5),  th(8),  th(1),  th(4),  th(7),      0,      0,      0;
        th(3),  th(6),  th(9),      0,      0,      0,  th(1),  th(4),  th(7);
            0,      0,      0,  th(3),  th(6),  th(9),  th(2),  th(5),  th(8)];
    end
end
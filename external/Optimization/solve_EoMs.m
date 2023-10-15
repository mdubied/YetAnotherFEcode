% solve_EoMs
%
% Synthax:
% TI_NL_PROM = solve_EoMs(V,ROM_Assembly,tailProperties,spineProperties,actuTop,actuBottom,h,tmax)
%
% Description: Computes the solutions for eta,dot{eta}, ddot{eta} for 
% [0,tmax] for a given ROM assembly and time step
%
% INPUTS: 
% (1) V:                    ROB   
% (2) ROM_Assembly:         ROM assembly           
% (3) tailProperties:       properties of the tail pressure force
%                           (matrices, tail elements etc.)
% (4) spineProperties:      properties of the spine change in momentum
%                           (tensor, spine elements etc.)
% (5) dragProperties:       properties of the drag force
% (6) actuTop:              vectors and matrices related to the actuation
%                           muscle at the top
% (7) actuBottom:           vectors and matrices related to the actuation
%                           muscle at the bottom
% (8) h:                    time step for time integration
% (9) tmax:                 simulation for [0,tmax]
%
% OUTPUTS:
% (1) TI_NL_PROM:           struct containing the solutions and related
%                           information
%     
% Last modified: 15/10/2023, Mathieu Dubied, ETH Zurich

function TI_NL_ROM = solve_EoMs(V,ROM_Assembly,tailProperties,spineProperties,actuTop,actuBottom,tmax)

    % SIMULATION PARAMETERS AND ICs _______________________________________
    eta0 = zeros(size(V,2),1);
    etad0 = zeros(size(V,2),1);
    etadd0 = zeros(size(V,2),1);

    % FORCES: ACTUATION, REACTIVE FORCE, DRAG FORCE _______________________
    % actuation properties
    B1T = actuTop.B1;
    B1B = actuBottom.B1;
    B2T = actuTop.B2;
    B2B = actuBottom.B2;
    k=1;

    % tail pressure force properties
    A = tailProperties.A;
    B = tailProperties.B;
    R = [0 -1 0 0 0 0;
         1 0 0 0 0 0;
         0 0 0 -1 0 0;
         0 0 1 0 0 0;
         0 0 0 0 0 1;
         0 0 0 0 -1 0];     % 90 degrees rotation counterclock-wise
    wTail = tailProperties.w;
    VTail = tailProperties.V;
    mTilde = tailProperties.mTilde;
    nodes = ROM_Assembly.Mesh.nodes;
    nodesInPos = V.'*reshape(nodes.',[],1);     % initial node position expressed in the ROM
    
    % external forces handles 
    actuSignalT = @(t) k/2*(1-(1+0.02*sin(t*2*pi)));    % to change below as well if needed
    actuSignalB = @(t) k/2*(1-(1-0.02*sin(t*2*pi)));
    
    fActu = @(t,q)  k/2*(1-(1+0.02*sin(t*2*pi)))*(B1T+B2T*q) + ...
                    k/2*(1-(1-0.02*sin(t*2*pi)))*(B1B+B2B*q);
                       
    fTail = @(q,qd)  0.5*mTilde*2*wTail^3*VTail.'*(dot(A*VTail*qd,R*B*VTail*(nodesInPos+q)).^2* ...
                        B*VTail*(nodesInPos+q));
    
    T = spineProperties.tensors.T;
    fSpine = @(q,qd,qdd) double(ttv(ttv(ttv(T,qdd,2),q+nodesInPos,2),q+nodesInPos,2) + ttv(ttv(ttv(T,qd,2),qd,2),q+nodesInPos,2) + ttv(ttv(ttv(T,qd,2),q+nodesInPos,2),qd,2));


    % NONLINEAR TIME INTEGRATION __________________________________________
    
    % instantiate object for nonlinear time integration
    TI_NL_ROM = ImplicitNewmark('timestep',h,'alpha',0.005,'MaxNRit',60,'RelTol',1e-6);
    
    % modal nonlinear Residual evaluation function handle
    Residual_NL_red = @(eta,etad,etadd,t)residual_reduced_nonlinear_actu_hydro(eta,etad,etadd, ...
        t,ROM_Assembly,fActu,fTail,fSpine,actuTop,actuBottom,actuSignalT,actuSignalB,tailProperties,spineProperties,R,nodesInPos);

    % time integration 
    TI_NL_ROM.Integrate(eta0,etad0,etadd0,tmax,Residual_NL_red);
    TI_NL_ROM.Solution.u = V * TI_NL_ROM.Solution.q; % get full order solution

end 

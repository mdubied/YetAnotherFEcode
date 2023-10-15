% solve_sensitivities
%
% Synthax:
% TI_sens = solve_sensitivities(V,xi_k,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,actuTop,actuBottom,etaSol,etadSol,etaddSol,h,tmax)
%
% Description: Computes the solutions for the sensitivity S for [0,tmax] 
% for a given PROM assembly, xi_k, and time step h
%
% INPUTS: 
% (1) V:                    ROB 
% (2) xi_k:                 current shape variation xi to consider
% (3) PROM_Assembly:        PROM assembly  
% (4) tensors_PROM          internal forces tensors
% (5) tailProperties:       properties of the tail pressure force
%                           (matrices, tail elements etc.)
% (6) spineProperties:      properties of the spine change in momentum
%                           (tensor, spine elements etc.)
% (7) dragProperties:       properties of the drag force
% (8) actuTop:              vectors and matrices related to the actuation
%                           muscle at the top
% (9) actuBottom:           vectors and matrices related to the actuation
%                           muscle at the bottom
% (10) etaSol:              solution for eta to consider
% (11) etadSol:             solution for dot{eta} to consider
% (12) etadSol:             solution for ddot{eta} to consider
% (13) h:                   time step for time integration
% (14) tmax:                simulation for [0,tmax]
%
% OUTPUTS:
% (1) TI_sens:              struct containing the solutions and related
%                           information
%     
%
% Last modified: 15/10/2023, Mathieu Dubied, ETH Zurich

function TI_sens = solve_sensitivities(V,xi_k,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,actuTop,actuBottom,etaSol,etadSol,etaddSol,h,tmax)        

    % SIMULATION PARAMETERS AND ICs _______________________________________
    s0 = zeros(size(V,2),size(xi_k,1));
    sd0 = zeros(size(V,2),size(xi_k,1));
    sdd0 = zeros(size(V,2),size(xi_k,1));
    
    % EVALUATE PARTIAL DERIVATIVES ALONG NOMINAL SOLUTION _________________
    
    % internal forces
    pd_fint = @(eta)DpROM_derivatives(eta,tensors_PROM);
    
    % actuation
    k=1;
    actuSignalT = @(t) k/2*(1-(1+0.02*sin(t*2*pi)));    % to change below as well if needed
    actuSignalB = @(t) k/2*(1-(1-0.02*sin(t*2*pi)));
    pd_actuTop = @(a)PROM_actu_derivatives(actuTop,a);
    pd_actuBottom = @(a)PROM_actu_derivatives(actuBottom,a);

    % tail pressure force
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
    UTail = tailProperties.U;
    mTilde = tailProperties.mTilde;
    nodes = PROM_Assembly.Mesh.nodes;
    nodesInPos = V.'*reshape(nodes.',[],1);     % initial node position expressed in the ROM
    xi=0;
    pd_tail = @(eta,etad) PROM_tail_pressure_derivatives(eta,etad,A,B,R,mTilde,wTail,nodesInPos,xi,VTail,UTail);

    % spine change in momentum
    spineTensors = spineProperties.tensors;
    pd_spine = @(eta,etad,etadd)PROM_spine_momentum_derivatives(eta,etad,etadd,nodesInPos,xi,spineTensors);
    
    % TIME INTEGRATION ____________________________________________________
    
    % instantiate object for time integration
    TI_sens = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true,'sens',true);
    
    % residual function handle
    Residual_sens = @(s,sd,sdd,t)residual_linear_sens(s,sd,sdd,t,PROM_Assembly, ...
        etaSol,etadSol,etaddSol, ...
        pd_fint,pd_tail,pd_spine,pd_actuTop,pd_actuBottom, ...
        actuSignalT,actuSignalB,h);

    % time integration
    TI_sens.Integrate(s0,sd0,sdd0,tmax,Residual_sens);
    
end
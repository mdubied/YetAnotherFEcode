% reduced_cost_function_w_constraints_TET4
%
% Synthax:
% Lr = reduced_cost_function_w_constraints_TET4(N,tailProperties,spineProperties,eta0,eta,etad,etadd,xiRebuild,xi,dr,AConstraint,bConstraint,barrierParam)
%
% Description:  Computes the cost function value in the ROB, considering
%               constraints on the shape variation parameters xi. Upper and
%               lower bounds inequality are considered. These inequality
%               constraints are included as (log) barrier functions. 
%
% INPUTS: 
% (1) N:                number of time steps             
% (2) tailProperties:   properties of the tail pressure force
%                       (matrices, tail elements etc.)
% (3) spineProperties:  properties of the spine change in momentum
%                       (tensor, spine elements etc.)
% (4) x0:               initial node position in FOM
% (5) eta:              solution for the reduced state variables
% (6) etad:             solution for the reduced velocities
% (7) etadd:            solution for the reduced accelerations
% (8) xiRebuild:        current value for xi, after the last PROM rebuild
% (9) xi:               current value for xi, after first PROM build
% (10) dr:              reduced forward swimming direction vector
% (11)-(12) A, b:       constraints on xi of the form Axi<b  
% (12) barrierParam:    parameter to scale (1/barrierParam) the barrier functions      
%                   
%
% OUTPUTS:
% (1) Lr:   reduced cost function value    
%     
%
% Last modified: 24/10/2023, Mathieu Dubied, ETH Zurich

function Lr = reduced_cost_function_w_constraints_TET4(N,tailProperties,spineProperties,x0,eta,etad,etadd,xiRebuild,xi,dr,AConstraint,bConstraint,barrierParam,V)
    Lr = 0;
    nConstraints = size(bConstraint);

    % TAIL PRESSURE FORCE _________________________________________________ 
    A = tailProperties.A;
    B = tailProperties.B;
    R = tailProperties.R;     
    wTail = tailProperties.w;
    VTail = tailProperties.V;
    UTail = tailProperties.U;
    z0 = tailProperties.z;
    iDOFs = tailProperties.iDOFs;
    x0Tail = x0(iDOFs);
    Uz = tailProperties.Uz;

    xDir = zeros(size(V,1),1);
    xDir(1:3:end) = 1;

    mTilde = 0.25*pi*1000*((z0+Uz*xiRebuild)*2)^2;
                
    fTail = @(q,qd)  0.5*mTilde*2*wTail^3*VTail.'*(dot(A*VTail*qd,R*B*(x0Tail+UTail*xi+VTail*q))).^2* ...
                        B*(x0Tail+UTail*xi+VTail*q);
    
    % SPINE CHANGE IN MOMENTUM ____________________________________________
    Txx = spineProperties.tensors.Txx.f0;
    TxV = spineProperties.tensors.TxV.f0;
    TVx = spineProperties.tensors.TVx.f0;
    TVV = spineProperties.tensors.TVV.f0;
    
    fSpine = @(q,qd,qdd) double(Txx)*qdd ...
        + double(ttv(ttv(TxV,q,3),qdd,2) ...
        + ttv(ttv(TVx,q,3),qdd,2) ...
        + ttv(ttv(ttv(TVV,q,4),q,3),qdd,2) ...
        + ttv(ttv(TVx,qd,3),qd,2) ...
        + ttv(ttv(ttv(TVV,q,4),qd,3),qd,2) ...
        + ttv(ttv(TxV,qd,3),qd,2) ...
        + ttv(ttv(ttv(TVV,qd,4),q,3),qd,2));
   
    LwoB = 0;
    
    for t=1:N-2
        eta_i = eta(:,t);
        etad_i = etad(:,t);
        etadd_i = etadd(:,t);

        % evaluate total hydrodynamic force
        % fhydro = fTail(eta_i,etad_i) + fSpine(eta_i,etad_i,etadd_i);
        
        % constraints (log barriers) to be included in the cost function
        logBarrierInTimeStep = 0;
        if nConstraints ~= 0
            for i = 1:nConstraints 
                logBarrierInTimeStep = logBarrierInTimeStep - 1/barrierParam*log(-AConstraint(i,:)*xi+bConstraint(i));
            end
        end

        % final cost function at time step t
        % Lr = Lr - dr'*fhydro + logBarrierInTimeStep;
        %Lr = [Lr, -dr'*fhydro + logBarrierInTimeStep];
        % LwoB = LwoB -dr'*fhydro;

        
        Lr = Lr - xDir.'*V*eta_i + logBarrierInTimeStep;
        LwoB = LwoB - xDir.'*V*eta_i;
    end

    % print cost function without part stemming from barrier functions
    fprintf('Drag minimization: %.4f\n',LwoB)
    fprintf('Full cost (with barrier): %.4f\n',Lr)
end
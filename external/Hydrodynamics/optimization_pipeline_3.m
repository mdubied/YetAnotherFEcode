% optimization_pipeline_3
%
% Synthax:
% [xiStar,LrEvo] = optimization_pipeline_3(MeshNominal,nodes,elements,U,d,h,tmax,FORMULATION,VOLUME,USEJULIA)
%
% Description: Implementation of the optimization pipeline 3 presented in
% the paper, based on Newton's method. 
%
% INPUTS: 
% (1) MeshNominal:  nominal mesh converted from Abaqus              
% (2) nodes:        nodes and their coordinates
% (3) elements:     elements and corresponding nodes
% (4) U:            shape variation basis
% (5) d:            forward swimming direction
% (6) h:            time step for time integration
% (7) tmax:         simulation for [0,tmax]
% (8) FORMULATION:  order of the Neumann approximation (N0/N1/N1t)
% (9) VOLUME:       integration over defected (1) or nominal volume (0)
% (10) USEJULIA:    use of JULIA (1) for the computation of internal forces
%                   tensors
%
% OUTPUTS:
% (1) xi_star:      optimal shape parameter(s) (scalar or vector)
% (2) LrEvo:        evolution of the cost function values
%     
%
% Additional notes: The functions "reduce_vector" and "cost_function" used
% by optimization_pipeline_3 are defined below.
%
%
% Last modified: 17/12/2022, Mathieu Dubied, ETH Zürich

function [xiStar,LrEvo] = optimization_pipeline_3(MeshNominal,nodes,elements,U,d,h,tmax,FORMULATION,VOLUME,USEJULIA)
    
    % STEP 1: set xi_0 = 0 ________________________________________________
    fprintf('____________________\n')
    fprintf('STEP 1\n')
    f = waitbar(0,'Step 1 ...','Name','Optimization pipeline P3');

    xi_k = 0;

    % STEP 2: mesh the structure and build a PROM _________________________
    fprintf('____________________\n')
    fprintf('STEP 2\n')
    waitbar(.05,f,'Step 2 ...');

    [V,PROM_Assembly,tensors_PROM,tensors_hydro_PROM] = ...
        build_PROM(MeshNominal,nodes,elements,U,FORMULATION,VOLUME,USEJULIA);
        
    % STEP 3: solve EoMs to get nominal solution eta and dot{eta} _________
    fprintf('____________________\n')
    fprintf('STEP 3\n')
    waitbar(.15,f,'Step 3 ...');
    
    TI_NL_PROM = solve_EoMs(V,PROM_Assembly,tensors_hydro_PROM,h,tmax);

    % STEP 4: solve sensitivity equation to get S _________________________
    fprintf('____________________\n')
    fprintf('STEP 4\n') 
    waitbar(.60,f,'Step 4 ...');

    TI_sens = solve_sensitivities(V,xi_k,PROM_Assembly,tensors_PROM, ...
        tensors_hydro_PROM,TI_NL_PROM.Solution.q,TI_NL_PROM.Solution.qd, ...
        TI_NL_PROM.Solution.qdd,h,tmax);
    
    % STEP 5-13: OPTIMIZATION LOOP ________________________________________
    fprintf('____________________\n')
    fprintf('STEP 5-13\n') 
    waitbar(.70,f,'Step 5-13 ...');
    % 
    close(f)

    eta = TI_NL_PROM.Solution.q;
    etad = TI_NL_PROM.Solution.qd;
    S = TI_sens.Solution.q;
    Sd = TI_sens.Solution.qd;
    eta_k = eta;
    etad_k = etad;

    N = size(eta,2);
    dr = reduce_vector(d,V);
    Lr = cost_function(N,tensors_hydro_PROM,eta,etad,dr);
    LrEvo = Lr;

    for k = 1:2
        fprintf('Optimization loop iteration: k= %d\n',k-1)
        % step 7
        eta_k = eta_k + S*xi_k; 
        % step 8
        nablaLr = gradient_cost_function(dr,xi_k,eta_k,etad_k,S,Sd,tensors_hydro_PROM);
        LrEvo = [LrEvo, cost_function(N,tensors_hydro_PROM,eta_k,etad_k,dr)];
        % step 9 and 10
        xi_k = xi_k - 0.8*nablaLr;
    end
    xiStar = xi_k;

end


% reduce_vector
% 
% Description: Helper function to compute the reduce version of the vector
% d, i.e., dr
%
function dr = reduce_vector(d,V)
    n = size(V,1);
    m = size(V,2);
    dr = zeros(m,1);
    for i=1:6:n
        Ve = V(i:i+5,:);
        de = [d;d;d];
        dr = dr + Ve.'*de;  
    end
end

% cost_function
%
% Description: Computes the cost function for a give eta, etad,dr and PROM
%
function Lr = cost_function(N,tensors_hydro_PROM,eta,etad,dr)
    Lr = 0;
    for i=1:N-2
        eta_i = eta(:,i);
        etad_i = etad(:,i);
        fhydro = double(tensors_hydro_PROM.Tr1) + ...
        double(tensors_hydro_PROM.Tru2*eta_i) + double(tensors_hydro_PROM.Trudot2*etad_i) + ...
        double(ttv(ttv(tensors_hydro_PROM.Truu3,eta_i,3), eta_i,2)) + ...
        double(ttv(ttv(tensors_hydro_PROM.Truudot3,etad_i,3), eta_i,2)) + ...
        double(ttv(ttv(tensors_hydro_PROM.Trudotudot3,etad_i,3), etad_i,2));
       
        Lr = Lr - dr'*fhydro;
    end
end

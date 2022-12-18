% optimization_pipeline_1
%
% Synthax:
% [xiStar,LrEvo] = optimization_pipeline_1(MeshNominal,nodes,elements,U,d,h,tmax,FORMULATION,VOLUME,USEJULIA)
%
% Description: Implementation of the optimization pipeline 1 presented in
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
% Additional notes: none
%
%
% Last modified: 17/12/2022, Mathieu Dubied, ETH Zürich
function [xiStar,LrEvo] = optimization_pipeline_1(myElementConstructor,nset,nodes,elements,U,d,h,tmax,FORMULATION,VOLUME,USEJULIA)
    
    % STEP 1: set xi_0 = 0 ________________________________________________
    fprintf('____________________\n')
    fprintf('STEP 1\n')
    %f = waitbar(0,'Step 1 ...','Name','Optimization pipeline P1');

    xi_k = 0;
    
    % STEP 2: set S_0 = 0 _________________________________________________
    fprintf('____________________\n')
    fprintf('STEP 2\n')
    %waitbar(0.02,'Step 2 ...');

    S_0 = 0;

    % STEP 3-12: optimization loop ________________________________________
    fprintf('____________________\n')
    fprintf('STEP 3-12\n')
    %waitbar(0.5,'Step 2 ...');


    for k = 1:10
        fprintf('Optimization loop iteration: k = %d\n',k-1)
        % STEP 5: mesh the structure and build a PROM _____________________
        % update defected mesh nodes
        df = U*xi_k;                       % displacement fields introduced by defects
        ddf = [df(1:2:end) df(2:2:end)]; 
        nodes_defected = nodes + ddf;    % nominal + d ---> defected 
        DefectedMesh = Mesh(nodes_defected);
        DefectedMesh.create_elements_table(elements,myElementConstructor);
        DefectedMesh.set_essential_boundary_condition([nset{1} nset{2}],1:2,0)

        [V,PROM_Assembly,tensors_PROM,tensors_hydro_PROM] = ...
        build_PROM(DefectedMesh,nodes_defected,elements,U,FORMULATION,VOLUME,USEJULIA);
        
        dr = reduced_constant_vector(d,V);

        % STEP 6: solve EoMs to get nominal solution eta and dot{eta} _____
        TI_NL_PROM = solve_EoMs(V,PROM_Assembly,tensors_hydro_PROM,h,tmax);

        eta_k = TI_NL_PROM.Solution.q;
        etad_k = TI_NL_PROM.Solution.qd;

        N = size(eta_k,2);
        
        % STEP 7: solve sensitivity equation to get S _____________________
        TI_sens = solve_sensitivities(V,xi_k,PROM_Assembly,tensors_PROM, ...
        tensors_hydro_PROM,TI_NL_PROM.Solution.q,TI_NL_PROM.Solution.qd, ...
        TI_NL_PROM.Solution.qdd,h,tmax);

        S_k = TI_sens.Solution.q;
        Sd_k = TI_sens.Solution.qd;

        % STEP 8: evaluate gradient _______________________________________
        nablaLr = gradient_cost_function(dr,xi_k,eta_k,etad_k,S_k,Sd_k,tensors_hydro_PROM);

        % STEP 9-10: update xi_k __________________________________________
        if k==1
            LrEvo = reduced_cost_function(N,tensors_hydro_PROM,eta_k,etad_k,dr);
        else
            LrEvo = [LrEvo, reduced_cost_function(N,tensors_hydro_PROM,eta_k,etad_k,dr)];
        end
        xi_k = xi_k - 0.8*nablaLr;
    end
    xiStar = xi_k;

end



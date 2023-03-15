% optimization_pipeline_1
%
% Synthax:
% [xiStar,LrEvo] = optimization_pipeline_1(myElementConstructor,nset,nodes,elements,U,d,h,tmax,A,b,varargin)
%
% Description: Implementation of the optimization pipeline 1 presented in
% the paper, based on Newton's method. 
%
% INPUTS: 
% (1) myElementConstructor: 
% (2) nset:         boundary conditions
% (3) nodes:        set of nodes of nominal (initial) mesh
% (4) elements      set of elements of nominal (initial) mesh
% (5) U:            shape variation basis
% (6) d:            forward swimming direction
% (7) h:            time step for time integration
% (8) tmax:         simulation for [0,tmax]
% (9)-(11)          constraints on xi of the form Axi<b
%
% possible additional name-value pair arguments
% (12) maxIteration:maximum number of iterations
% (13) convCrit:    convergence criterium. Norm between two successive
%                   optimal paramter vectors
% (14) FORMULATION: order of the Neumann approximation (N0/N1/N1t)
% (15) VOLUME:      integration over defected (1) or nominal volume (0)
% (16) USEJULIA:    use of JULIA (1) for the computation of internal forces
%                   tensors
%
% OUTPUTS:
% (1) xi_star:      optimal shape parameter(s) (scalar or vector)
% (2) LrEvo:        evolution of the cost function values
%     
%
% Last modified: 15/03/2023, Mathieu Dubied, ETH Zurich
function [xiStar,xiEvo,LrEvo] = optimization_pipeline_1(myElementConstructor,nset,nodes,elements,U,d,h,tmax,A,b,varargin)
    
    % parse input
    [maxIteration,convCrit,FORMULATION,VOLUME,USEJULIA,FOURTHORDER] = parse_inputs(varargin{:});

    % STEP 1: set xi_0 = 0 ________________________________________________
    fprintf('____________________\n')
    fprintf('STEP 1\n')
    %f = waitbar(0,'Step 1 ...','Name','Optimization pipeline P1');

    xi_k = 0;
    xiEvo = xi_k;
    

    % STEP 2-12: optimization loop ________________________________________
    fprintf('____________________\n')
    fprintf('STEP 2-11\n')
    %waitbar(0.5,'Step 2 ...');


    for k = 1:maxIteration
        fprintf('Optimization loop iteration: k = %d\n',k-1)
        % STEP 4: mesh the structure and build a PROM _____________________
        % update defected mesh nodes
        df = U*xi_k;                       % displacement fields introduced by defects
        ddf = [df(1:2:end) df(2:2:end)]; 
        nodes_defected = nodes + ddf;    % nominal + d ---> defected 
        DefectedMesh = Mesh(nodes_defected);
        DefectedMesh.create_elements_table(elements,myElementConstructor);
        DefectedMesh.set_essential_boundary_condition([nset{1} nset{2}],1:2,0)

        [V,PROM_Assembly,tensors_PROM,tensors_hydro_PROM] = ...
        build_PROM(DefectedMesh,nodes_defected,elements,U,FORMULATION,VOLUME,USEJULIA,FOURTHORDER);
        
        dr = reduced_constant_vector(d,V);

        % STEP 5: solve EoMs to get nominal solution eta and dot{eta} _____
        TI_NL_PROM = solve_EoMs(V,PROM_Assembly,tensors_hydro_PROM,h,tmax);

        eta_k = TI_NL_PROM.Solution.q;
        etad_k = TI_NL_PROM.Solution.qd;

        N = size(eta_k,2);
        
        % STEP 6: solve sensitivity equation to get S _____________________
        TI_sens = solve_sensitivities(V,xi_k,PROM_Assembly,tensors_PROM, ...
        tensors_hydro_PROM,TI_NL_PROM.Solution.q,TI_NL_PROM.Solution.qd, ...
        TI_NL_PROM.Solution.qdd,h,tmax,FOURTHORDER);

        S_k = TI_sens.Solution.q;
        Sd_k = TI_sens.Solution.qd;

        % STEP 7: evaluate gradient _______________________________________
        nablaLr = gradient_cost_function_w_constraints(dr,xi_k,eta_k,etad_k,S_k,Sd_k,A,b,tensors_hydro_PROM,FOURTHORDER);

        % STEP 8-9: update xi_k __________________________________________
        if k==1
            LrEvo = reduced_cost_function_w_constraints(N,tensors_hydro_PROM,eta_k,etad_k,xi_k,dr,A,b);
        else 
            LrEvo = [LrEvo, reduced_cost_function_w_constraints(N,tensors_hydro_PROM,eta_k,etad_k,xi_k,dr,A,b)];
        end
        xi_k = xi_k - 0.8*nablaLr;
        xiEvo = [xiEvo,xi_k];
        
        if norm(xiEvo(end)-xiEvo(end-1))<0.001
            fprintf('Convergence criterion of %.2g fulfilled\n',convCrit)
            break
        end
        if k == maxIteration
            fprintf('Maximum number of %d iterations reached\n',maxIteration)
        end
    end

    xiStar = xi_k;

end

% parse input
function [maxIteration,convCrit,FORMULATION,VOLUME,USEJULIA,FOURTHORDER] = parse_inputs(varargin)
defaultMaxIteration = 10;
defaultConvCrit = 0.001;
defaultFORMULATION = 'N1';
defaultVOLUME = 1;
defaultUSEJULIA = 0; 
defaultFOURTHORDER = 0;

p = inputParser;
addParameter(p,'maxIteration',defaultMaxIteration, @(x)validateattributes(x, ...
                {'numeric'},{'nonempty','integer','positive'}) );
addParameter(p,'convCrit',defaultConvCrit,@(x)validateattributes(x, ...
                {'numeric'},{'nonempty','positive'}) );
addParameter(p,'FORMULATION',defaultFORMULATION,@(x)validateattributes(x, ...
                {'char'},{'nonempty'}))
addParameter(p,'VOLUME',defaultVOLUME,@(x)validateattributes(x, ...
                {'numeric'},{'nonempty'}) );
addParameter(p,'USEJULIA',defaultUSEJULIA,@(x)validateattributes(x, ...
                {'numeric'},{'nonempty'}) );
addParameter(p,'FOURTHORDER',defaultFOURTHORDER,@(x)validateattributes(x, ...
                {'numeric'},{'nonempty','integer','positive'}) );

parse(p,varargin{:});

maxIteration = p.Results.maxIteration;
convCrit = p.Results.convCrit;
FORMULATION = p.Results.FORMULATION;
VOLUME = p.Results.VOLUME;
USEJULIA = p.Results.USEJULIA;
FOURTHORDER = p.Results.FOURTHORDER;

end

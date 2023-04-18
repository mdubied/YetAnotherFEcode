% optimization_pipeline_2
%
% Synthax:
% [xiStar,LrEvo] = optimization_pipeline_2(MeshNominal,nodes,elements,U,d,h,tmax,A,b,varargin)
%
% Description: Implementation of the optimization pipeline 2 presented in
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
% (9)-(10)          constraints on xi of the form Axi<b
%
% possible additional name-value pair arguments
% (11) maxIteration:maximum number of iterations
% (12) convCrit:    convergence criterium. Norm between two successive
%                   optimal paramter vectors
% (13) FORMULATION: order of the Neumann approximation (N0/N1/N1t)
% (14) VOLUME:      integration over defected (1) or nominal volume (0)
% (15) USEJULIA:    use of JULIA (1) for the computation of internal forces
%                   tensors
% (16) barrierParam:parameter to scale the barrier function for the 
%                   constraints (1/barrierParam)
% (17) gStepSize:   step size used in the gradient descent algorithm
%
%
% OUTPUTS:
% (1) xi_star:      optimal shape parameter(s) (scalar or vector)
% (2) LrEvo:        evolution of the cost function values
%     
%
% Last modified: 18/04/2023, Mathieu Dubied, ETH Zurich

function [xiStar,xiEvo,LrEvo] = optimization_pipeline_2(MeshNominal,nodes,elements,U,d,h,tmax,A,b,varargin)
    
     % parse input
    [maxIteration,convCrit,barrierParam,gStepSize,FORMULATION,VOLUME,USEJULIA,FOURTHORDER,ACTUATION] = parse_inputs(varargin{:});
    
    % STEP 1: set xi_0 = 0 ________________________________________________
    fprintf('____________________\n')
    fprintf('STEP 1\n')

    xi_k = zeros(size(U,2),1);
    xiEvo = xi_k;
    
    % STEP 2: mesh the structure and build a PROM _________________________
    fprintf('____________________\n')
    fprintf('STEP 2\n')    
    
    [V,PROM_Assembly,tensors_PROM,tensors_hydro_PROM, tensors_topMuscle_PROM, tensors_bottomMuscle_PROM] = ...
        build_PROM(MeshNominal,nodes,elements,U,FORMULATION,VOLUME,USEJULIA,FOURTHORDER,ACTUATION);

    
    % STEP 3: solve EoMs to get nominal solution eta_0 and dot{eta}_0 _____
    fprintf('____________________\n')
    fprintf('STEP 2\n')

    if ACTUATION
        TI_NL_PROM = solve_EoMs(V,PROM_Assembly,tensors_hydro_PROM,h,tmax,...
            'ACTUATION', ACTUATION,'topMuscle',tensors_topMuscle_PROM,'bottomMuscle',tensors_bottomMuscle_PROM);
    else
        TI_NL_PROM = solve_EoMs(V,PROM_Assembly,tensors_hydro_PROM,h,tmax);
    end

    eta = TI_NL_PROM.Solution.q;
    etad = TI_NL_PROM.Solution.qd;
      
    eta_k = eta;
    etad_k = etad;
    S_k = zeros(size(eta,1),size(xi_k,1));

    % STEP 4-12: optimization loop ________________________________________
    fprintf('____________________\n')
    fprintf('STEP 3-12\n')

    N = size(eta,2);
    dr = reduced_constant_vector(d,V);

    for k = 1:maxIteration
        fprintf('Optimization loop iteration: k = %d\n',k-1)
        % STEP 6: approximate eta_k _____________________
        eta_k = eta_k + S_k*xi_k; 

        % STEP 7: solve sensitivity equation to get S _____________________
        if ACTUATION
            TI_sens = solve_sensitivities(V,xi_k,PROM_Assembly,tensors_PROM, ...
                tensors_hydro_PROM,TI_NL_PROM.Solution.q,TI_NL_PROM.Solution.qd, ...
                TI_NL_PROM.Solution.qdd,h,tmax,FOURTHORDER,...
                    'ACTUATION', ACTUATION,'topMuscle',tensors_topMuscle_PROM,'bottomMuscle',tensors_bottomMuscle_PROM);
        else
             TI_sens = solve_sensitivities(V,xi_k,PROM_Assembly,tensors_PROM, ...
                tensors_hydro_PROM,TI_NL_PROM.Solution.q,TI_NL_PROM.Solution.qd, ...
                TI_NL_PROM.Solution.qdd,h,tmax,FOURTHORDER);
        end


        S_k = TI_sens.Solution.q;
        Sd_k = TI_sens.Solution.qd;

        % STEP 8: evaluate gradient _______________________________________
        nablaLr = gradient_cost_function_w_constraints(dr,xi_k,xi_k,eta_k,etad_k,S_k,Sd_k,A,b,barrierParam,tensors_hydro_PROM,FOURTHORDER);

        % STEP 9-10: update xi_k __________________________________________
        if k==1
            LrEvo = reduced_cost_function_w_constraints(N,tensors_hydro_PROM,eta_k,etad_k,xi_k,xi_k,dr,A,b,barrierParam);
        else
            LrEvo = [LrEvo, reduced_cost_function_w_constraints(N,tensors_hydro_PROM,eta_k,etad_k,xi_k,xi_k,dr,A,b,barrierParam)];
        end
        xi_k = xi_k - gStepSize*nablaLr
        xiEvo = [xiEvo,xi_k];
        
        % possible exit conditions
        if size(xi_k,1) >1
            if norm(xiEvo(:,end)-xiEvo(:,end-1))<convCrit
                fprintf('Convergence criterion of %.2g fulfilled\n',convCrit)
                break
            end
        else
            if norm(xiEvo(end)-xiEvo(end-1))<convCrit
                fprintf('Convergence criterion of %.2g fulfilled\n',convCrit)
                break
            end
        end

        if k == maxIteration
            fprintf('Maximum number of %d iterations reached\n',maxIteration)
        end
    end

    xiStar = xi_k;

end

% parse input
function [maxIteration,convCrit,barrierParam,gStepSize,FORMULATION,VOLUME,USEJULIA,FOURTHORDER,ACTUATION] = parse_inputs(varargin)
defaultMaxIteration = 50;
defaultConvCrit = 0.001;
defaultBarrierParam = 500;
defaultGStepSize = 0.1;
defaultFORMULATION = 'N1';
defaultVOLUME = 1;
defaultUSEJULIA = 0; 
defaultFOURTHORDER = 0;
defaultACTUATION = 0;
p = inputParser;
addParameter(p,'maxIteration',defaultMaxIteration, @(x)validateattributes(x, ...
                {'numeric'},{'nonempty','integer','positive'}) );
addParameter(p,'convCrit',defaultConvCrit,@(x)validateattributes(x, ...
                {'numeric'},{'nonempty','positive'}) );
addParameter(p,'barrierParam',defaultBarrierParam,@(x)validateattributes(x, ...
                {'numeric'},{'nonempty','positive'}) );
addParameter(p,'gStepSize',defaultGStepSize,@(x)validateattributes(x, ...
                {'numeric'},{'nonempty','positive'}) );
addParameter(p,'FORMULATION',defaultFORMULATION,@(x)validateattributes(x, ...
                {'char'},{'nonempty'}))
addParameter(p,'VOLUME',defaultVOLUME,@(x)validateattributes(x, ...
                {'numeric'},{'nonempty'}) );
addParameter(p,'USEJULIA',defaultUSEJULIA,@(x)validateattributes(x, ...
                {'numeric'},{'nonempty'}) );
addParameter(p,'FOURTHORDER',defaultFOURTHORDER,@(x)validateattributes(x, ...
                {'numeric'},{'nonempty','integer','positive'}) );
addParameter(p,'ACTUATION',defaultACTUATION,@(x)validateattributes(x, ...
                {'numeric'},{'nonempty','integer','positive'}) );

parse(p,varargin{:});

maxIteration = p.Results.maxIteration;
convCrit = p.Results.convCrit;
barrierParam = p.Results.barrierParam;
gStepSize = p.Results.gStepSize;
FORMULATION = p.Results.FORMULATION;
VOLUME = p.Results.VOLUME;
USEJULIA = p.Results.USEJULIA;
FOURTHORDER = p.Results.FOURTHORDER;
ACTUATION = p.Results.ACTUATION;

end
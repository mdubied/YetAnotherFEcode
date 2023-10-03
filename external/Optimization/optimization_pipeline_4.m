% optimization_pipeline_4
%
% Synthax:
% [xiStar,LrEvo] = optimization_pipeline_4(MeshNominal,nodes,elements,U,d,h,tmax,A,b,varargin)
%
% Description: Implementation of the optimization pipeline 4 presented in
% the paper, based on Newton's method. 
%
% INPUTS: 
% (1) myElementConstructor: 
% (2) nset:         boundary conditions             
% (3) nodes:        nodes and their coordinates
% (4) elements:     elements and corresponding nodes
% (5) U:            shape variation basis
% (6) d:            forward swimming direction
% (7) h:            time step for time integration
% (8) tmax:         simulation for [0,tmax]
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
% (16) ACTUATION:   describes cases with actuation. 0 is w/o actuation, 1
%                   with actuation for the 2D case discussed in the paper
% (17) barrierParam:parameter to scale the barrier function for the 
%                   constraints (1/barrierParam)
% (18) gStepSize:   step size used in the gradient descent algorithm
% (19) nRebuild:    number of step between each re-build of a PROM
%
% OUTPUTS:
% (1) xi_star:      optimal shape parameter(s) (scalar or vector)
% (2) LrEvo:        evolution of the cost function values
%     
%
% Last modified: 02/04/2023, Mathieu Dubied, ETH Zurich

function [xiStar,xiEvo,LrEvo] = optimization_pipeline_4(myElementConstructor,nset,nodes,elements,U,d,h,tmax,A,b,varargin)

    % parse input
    [maxIteration,convCrit,barrierParam,gStepSize,nRebuild,FORMULATION,VOLUME,USEJULIA,FOURTHORDER,ACTUATION] = parse_inputs(varargin{:});
    
    % STEP 1: set xi_0 = 0 ________________________________________________
    fprintf('____________________\n')
    fprintf('STEP 1\n')

    xi_k = zeros(size(U,2),1);
    xiRebuild_k = zeros(size(U,2),1);
    xiEvo = xi_k;

    % STEP 2: mesh the structure and build a PROM _________________________
    fprintf('____________________\n')
    fprintf('STEP 2\n')

    % create nominal mesh
    MeshNominal = Mesh(nodes);
    MeshNominal.create_elements_table(elements,myElementConstructor);
    MeshNominal.set_essential_boundary_condition([nset{1} nset{2}],1:2,0)

    [V,PROM_Assembly,tensors_PROM,tensors_hydro_PROM, tensors_topMuscle_PROM, tensors_bottomMuscle_PROM] = ...
        build_PROM(MeshNominal,nodes,elements,U,FORMULATION,VOLUME,USEJULIA,FOURTHORDER,ACTUATION);

    
    % STEP 3: solve EoMs to get nominal solution eta and dot{eta} _________
    fprintf('____________________\n')
    fprintf('STEP 3\n')
    
    if ACTUATION
        TI_NL_PROM = solve_EoMs(V,PROM_Assembly,tensors_hydro_PROM,h,tmax,...
            'ACTUATION', ACTUATION,'topMuscle',tensors_topMuscle_PROM,'bottomMuscle',tensors_bottomMuscle_PROM);
    else
        TI_NL_PROM = solve_EoMs(V,PROM_Assembly,tensors_hydro_PROM,h,tmax);
    end
    

    % STEP 4: solve sensitivity equation to get S _________________________
    fprintf('____________________\n')
    fprintf('STEP 4\n') 
    
    if ACTUATION
        TI_sens = solve_sensitivities(V,xiRebuild_k,PROM_Assembly,tensors_PROM, ...
            tensors_hydro_PROM,TI_NL_PROM.Solution.q,TI_NL_PROM.Solution.qd, ...
            TI_NL_PROM.Solution.qdd,h,tmax,FOURTHORDER,...
                'ACTUATION', ACTUATION,'topMuscle',tensors_topMuscle_PROM,'bottomMuscle',tensors_bottomMuscle_PROM);
    else
         TI_sens = solve_sensitivities(V,xiRebuild_k,PROM_Assembly,tensors_PROM, ...
            tensors_hydro_PROM,TI_NL_PROM.Solution.q,TI_NL_PROM.Solution.qd, ...
            TI_NL_PROM.Solution.qdd,h,tmax,FOURTHORDER);
    end

    
    % STEP 5-13: OPTIMIZATION LOOP ________________________________________
    fprintf('____________________\n')
    fprintf('STEP 5-13\n') 

    eta = TI_NL_PROM.Solution.q;
    etad = TI_NL_PROM.Solution.qd;
    S = TI_sens.Solution.q;
    Sd = TI_sens.Solution.qd;
    eta_k = eta;
    etad_k = etad;

    N = size(eta,2);
    dr = reduced_constant_vector(d,V);
    Lr = reduced_cost_function_w_constraints(N,tensors_hydro_PROM,eta,etad,xiRebuild_k,xi_k,dr,A,b,barrierParam);
    LrEvo = Lr;

    for k = 1:maxIteration

        fprintf('Optimization loop iteration: k= %d\n',k-1)

        % possible rebuilding of a PROM
        if mod(k,nRebuild) == 0
             
            fprintf('PROM Rebuild ... \n')
            % update defected mesh nodes
            df = U*xi_k;                       % displacement fields introduced by defects
            ddf = [df(1:2:end) df(2:2:end)]; 
            nodes_defected = nodes + ddf;    % nominal + d ---> defected 
            svMesh = Mesh(nodes_defected);
            svMesh.create_elements_table(elements,myElementConstructor);
            svMesh.set_essential_boundary_condition([nset{1} nset{2}],1:2,0)
    
            [V,PROM_Assembly,tensors_PROM,tensors_hydro_PROM, tensors_topMuscle_PROM, tensors_bottomMuscle_PROM] = ...
            build_PROM(svMesh,nodes,elements,U,FORMULATION,VOLUME,USEJULIA,FOURTHORDER,ACTUATION);
            
            xiRebuild_k = zeros(size(U,2),1);   % reset local xi to 0 as we rebuild the ROM
            
            dr = reduced_constant_vector(d,V);
    
            % solve EoMs to get nominal solution eta and dot{eta}
            if ACTUATION
                TI_NL_PROM = solve_EoMs(V,PROM_Assembly,tensors_hydro_PROM,h,tmax,...
                    'ACTUATION', ACTUATION,'topMuscle',tensors_topMuscle_PROM,'bottomMuscle',tensors_bottomMuscle_PROM);
            else
                TI_NL_PROM = solve_EoMs(V,PROM_Assembly,tensors_hydro_PROM,h,tmax);
            end
    
            eta_k = TI_NL_PROM.Solution.q;
            etad_k = TI_NL_PROM.Solution.qd;
    
            N = size(eta_k,2);
            
            % solve sensitivity equation to get S
            if ACTUATION
                TI_sens = solve_sensitivities(V,xiRebuild_k,PROM_Assembly,tensors_PROM, ...
                    tensors_hydro_PROM,TI_NL_PROM.Solution.q,TI_NL_PROM.Solution.qd, ...
                    TI_NL_PROM.Solution.qdd,h,tmax,FOURTHORDER,...
                        'ACTUATION', ACTUATION,'topMuscle',tensors_topMuscle_PROM,'bottomMuscle',tensors_bottomMuscle_PROM);
            else
                 TI_sens = solve_sensitivities(V,xiRebuild_k,PROM_Assembly,tensors_PROM, ...
                    tensors_hydro_PROM,TI_NL_PROM.Solution.q,TI_NL_PROM.Solution.qd, ...
                    TI_NL_PROM.Solution.qdd,h,tmax,FOURTHORDER);
            end
    
    
            S = TI_sens.Solution.q;
            Sd = TI_sens.Solution.qd;
        else
            % step 7
            if size(xi_k,1)>1
                S=tensor(S);
                eta_k = eta_k + double(ttv(S,xiRebuild_k,2));
            else
                eta_k = eta_k + S*xiRebuild_k;
            end
        end 

        % step 8
        nablaLr = gradient_cost_function_w_constraints(dr,xiRebuild_k,xi_k,eta_k,etad_k,S,Sd,A,b,barrierParam,tensors_hydro_PROM,FOURTHORDER);
        LrEvo = [LrEvo, reduced_cost_function_w_constraints(N,tensors_hydro_PROM,eta_k,etad_k,xiRebuild_k,xi_k,dr,A,b,barrierParam)];
        % step 9 and 10
        xi_k = xi_k - gStepSize*nablaLr
        xiRebuild_k = xiRebuild_k - gStepSize*nablaLr;

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
function [maxIteration,convCrit,barrierParam,gStepSize,nRebuild,FORMULATION,VOLUME,USEJULIA,FOURTHORDER,ACTUATION] = parse_inputs(varargin)
defaultMaxIteration = 50;
defaultConvCrit = 0.001;
defaultBarrierParam = 500;
defaultGStepSize = 0.1;
defaultNRebuild = 10;
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
addParameter(p,'nRebuild',defaultNRebuild,@(x)validateattributes(x, ...
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
nRebuild = p.Results.nRebuild;
FORMULATION = p.Results.FORMULATION;
VOLUME = p.Results.VOLUME;
USEJULIA = p.Results.USEJULIA;
FOURTHORDER = p.Results.FOURTHORDER;
ACTUATION = p.Results.ACTUATION;

end
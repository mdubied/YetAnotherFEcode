% DO NOT USE
% optimization_pipeline_5
%
% Synthax:
% [xiStar,LrEvo] = optimization_pipeline_5(myElementConstructor,nset,nodes,elements,U,d,h,tmax,A,b,varargin)
%
% Description: Implementation of the optimization pipeline 5 presented in
% the paper, based on Newton's method. P5 is using a PFOM model.
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
% Additional note:  Cannot be used. As we need a PFOM model, the internal
%                   forces tensors need to be construction on the on FOM. 
%                   This will saturate the memory.
%
% Last modified: 02/04/2023, Mathieu Dubied, ETH Zurich
function [xiStar,xiEvo,LrEvo] = optimization_pipeline_5(myElementConstructor,nset,nodes,elements,U,d,h,tmax,A,b,varargin)
    
    % parse input
    [maxIteration,convCrit,barrierParam,gStepSize,FORMULATION,VOLUME,USEJULIA,FOURTHORDER,ACTUATION] = parse_inputs(varargin{:});

    if ACTUATION
        error('Optimization Pipeline P4 cannot handle actuation.')
    end

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
        nodes_shapeVaried = nodes + ddf;    % nominal + d ---> defected 
        svMesh = Mesh(nodes_shapeVaried);
        svMesh.create_elements_table(elements,myElementConstructor);
        svMesh.set_essential_boundary_condition([nset{1} nset{2}],1:2,0)


        % Assembly
        svAssembly = Assembly(svMesh);
        Msv = svAssembly.mass_matrix();
        nNodes = size(nodes,1);
        u0 = zeros( svMesh.nDOFs, 1);
        [Ksv,~] = svAssembly.tangent_stiffness_and_force(u0);
        % store matrices
        svAssembly.DATA.K = Ksv;
        svAssembly.DATA.M = Msv;

        % damping 
        alfa = 3.1;
        beta = 6.3*1e-6;
        Csv = alfa*Msv + beta*Ksv; % Rayleigh damping 
        svAssembly.DATA.C = Csv;

        % internal force as a function of xi
        FORMULATION = 'N1'; % N1/N1t/N0
        VOLUME = 1;         % integration over defected (1) or nominal volume (0)
        USEJULIA = 0;
        V=sparse(ones(svMesh.nDOFs));

        % will run out of memory!
        tensors_PFOM = reduced_tensors_DpROM(svAssembly, elements, ...
        V, U, FORMULATION, VOLUME, USEJULIA); %compute tensors
        
        % hydrodynamic tensors
        [~,~,skinElements, skinElementFaces] = getSkin2D(elements);
        vwater = [1;0.1];   % water velocity vector
        rho = 1;
        tensors_hydro_PFOM = unreduced_tensors_hydro_PFOM(svAssembly, elements, U, skinElements, skinElementFaces, vwater, rho);

%         [V,PROM_Assembly,tensors_PROM,tensors_hydro_PROM] = ...
%         build_PROM(svMesh,nodes_shapeVaried,elements,U,FORMULATION,VOLUME,USEJULIA,FOURTHORDER);
%         
        %dr = reduced_constant_vector(d,V);
        
        d = repmat(d',1,nNodes)';

        % STEP 5: solve EoMs to get nominal solution eta and dot{eta} _____
        TI_NL_PROM = solve_EoMs_FOM(svAssembly,tensors_hydro_PFOM,h,tmax);

        eta_k = TI_NL_PROM.Solution.q; % or u, unconstrained?
        etad_k = TI_NL_PROM.Solution.qd;

        N = size(eta_k,2);
        
        % STEP 6: solve sensitivity equation to get S _____________________
        TI_sens = solve_sensitivities(V,xi_k,svAssembly,tensors_PROM, ...
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
        xi_k = xi_k - gStepSize*nablaLr;
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

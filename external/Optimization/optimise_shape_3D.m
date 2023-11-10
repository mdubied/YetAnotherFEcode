% optimise_shape_3D
%
% Synthax:
% [xiStar,xiEvo,LrEvo] = optimise_shape_3D(myElementConstructor,nset,nodes,elements,U,d,h,tmax,A,b,varargin)
%
% Description: Implementation of the optimization pipeline 4 presented in
% the paper, based on Newton's method. 
%
% INPUTS: 
% (1) myElementConstructor: defines the type of element and material
%                           properties
% (2) nset:                 set of element to apply boundary conditions             
% (3) nodes:                nodes and their coordinates
% (4) elements:             elements described by their nodes
% (5) U:                    shape variation basis
% (6) d:                    forward swimming direction
% (7) h:                    time step for time integration
% (8) tmax:                 simulation for [0,tmax]
% (9)-(10) A,b              constraints on xi of the form Axi<b
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
% (18) nRebuild:    number of step between each re-build of a PROM
%
% OUTPUTS:
% (1) xiStar:       optimal shape parameter(s) (scalar or vector)
% (2) xiEvo:        evolution of the optimal shape parameter(s)
% (3) LrEvo:        evolution of the cost function values
%     
%
% Last modified: 30/10/2023, Mathieu Dubied, ETH Zurich

function [xiStar,xiEvo,LrEvo] = optimise_shape_3D(myElementConstructor,nset,nodes,elements,U,d,h,tmax,A,b,varargin)

    % parse input
    [maxIteration,convCrit,convCritCost,barrierParam,gStepSize,nRebuild,...
        rebuildThreshold,FORMULATION,VOLUME,USEJULIA] = parse_inputs(varargin{:});
    
    % NOMINAL SOLUTION ____________________________________________________
    fprintf('**************************************\n')
    fprintf('Solving for the nominal structure...\n')
    fprintf('**************************************\n')

    xi_k = zeros(size(U,2),1);
    xiRebuild_k = zeros(size(U,2),1);
    xiEvo = xi_k;

    % Mesh
            
    MeshNominal = Mesh(nodes);
    MeshNominal.create_elements_table(elements,myElementConstructor);
 
    for l=1:length(nset)
        MeshNominal.set_essential_boundary_condition([nset{l}],1:3,0)   
    end

    % build PROM
    fprintf('____________________\n')
    fprintf('Building PROM ... \n')

    mTilde = 10;
    [V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom] = ...
    build_PROM_3D(MeshNominal,nodes,elements,mTilde,U,USEJULIA,VOLUME,FORMULATION);      
    

    % Solve EoMs
    tic 
    fprintf('____________________\n')
    fprintf('Solving EoMs...\n') 
    % TI_NL_PROM = solve_EoMs(V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom,h,tmax);      
    TI_NL_PROM = solve_EoMs_and_sensitivities(V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom,h,tmax);                        
    toc

    uTail = zeros(3,tmax/h);
    timePlot = linspace(0,tmax-h,tmax/h);
    x0Tail = min(nodes(:,1));
    for a=1:tmax/h
        uTail(:,a) = V(tailProperties.tailNode*3-2:tailProperties.tailNode*3,:)*TI_NL_PROM.Solution.q(:,a);
    end

    figure
    subplot(2,1,1);
    plot(timePlot,x0Tail+uTail(1,:),'DisplayName','k=0')
    hold on
    xlabel('Time [s]')
    ylabel('x-position tail node')
    legend('Location','northwest')
    subplot(2,1,2);
    plot(timePlot,uTail(2,:),'DisplayName','k=0')
    hold on
    xlabel('Time [s]')
    ylabel('y-position tail node')
    legend('Location','southwest')
    drawnow


    % Solve sensitivity equation 
    % tic
    % fprintf('____________________\n')
    % fprintf('Solving sensitivity...\n') 
    % TI_sens = solve_sensitivities(V,xi_k,PROM_Assembly, ...
    %    tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom, ...
    %    TI_NL_PROM.Solution.q,TI_NL_PROM.Solution.qd,TI_NL_PROM.Solution.qdd, ...
    %    h,tmax);
    % toc
    
    % Retrieving solutions    
    eta = TI_NL_PROM.Solution.q;
    etad = TI_NL_PROM.Solution.qd;
    etadd = TI_NL_PROM.Solution.qdd;
    % S = TI_sens.Solution.q;
    % Sd = TI_sens.Solution.qd;
    % Sdd = TI_sens.Solution.qdd;

    S = TI_NL_PROM.Solution.s;
    Sd = TI_NL_PROM.Solution.sd;
    Sdd = TI_NL_PROM.Solution.sdd;
    eta_0k = TI_NL_PROM.Solution.q;
    eta_k = eta;
    etad_k = etad;
    etadd_k = etadd;

    % nodes = PROM_Assembly.Mesh.nodes;
    x0 = reshape(nodes.',[],1);
    
    % computing initial cost function value
    fprintf('____________________\n')
    fprintf('Computing cost function...\n') 
    N = size(eta,2);
    dr = reduced_constant_vector(d,V,3);
    Lr = reduced_cost_function_w_constraints_TET4(N,tailProperties,spineProperties,x0,eta,etad,etadd,xiRebuild_k,xi_k,dr,A,b,barrierParam,V);  
    LrEvo = Lr;
    nablaEvo = zeros(size(U,2),1);
    lastRebuild = 0;

    for k = 1:maxIteration
        fprintf('**************************************\n')
        fprintf('Optimization loop iteration: k= %d\n',k)
        fprintf('**************************************\n')

        % possible rebuilding of a PROM
        if check_cond_rebuild(k,lastRebuild,nRebuild,xiRebuild_k,rebuildThreshold,maxIteration)
            lastRebuild = k;

            % update defected mesh nodes
            df = U*xi_k;                       % displacement fields introduced by defects
            ddf = [df(1:3:end) df(2:3:end) df(3:3:end)]; 
            nodes_defected = nodes + ddf;    % nominal + d ---> defected 
            svMesh = Mesh(nodes_defected);
            svMesh.create_elements_table(elements,myElementConstructor);
            for l=1:length(nset)
                svMesh.set_essential_boundary_condition([nset{l}],1:3,0)   
            end

            % build PROM
            [V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom] = ...
                 build_PROM_3D(svMesh,nodes_defected,elements,mTilde,U,USEJULIA,VOLUME,FORMULATION);

            % Lx = abs(max(nodes(:,1))-min(nodes(:,1)));  % horizontal length of airfoil
            % Ly = abs(max(nodes(:,2))-min(nodes(:,2)));  % vertical length of airfoil
            % Lz = abs(max(nodes(:,3))-min(nodes(:,3)));  % vertical length of airfoil
            % f1 = figure('units','centimeters','position',[3 3 15 7],'name','Shape-varied mesh');
            % elementPlot = elements(:,1:4); hold on 
            % % PlotMesh(nodes, elementPlot, 0); 
            % v1 = reshape(U*xi_k, 3, []).';
            % S = 1;
            % hf=PlotFieldonDeformedMesh(nodes, elementPlot, v1, 'factor', S);
            % L = [Lx,Ly,Lz];
            % O = [-Lx,-Ly/2,-Lz/2];
            % plotcube(L,O,.05,[0 0 0]);
            % axis equal; grid on; box on; drawnow
            % set(f1,'PaperUnits','centimeters');
            % set(f1,'PaperPositionMode','auto');
            % % set(f1,'PaperSize',[10 3.5]); % Canvas Size
            % set(f1,'Units','centimeters');
                                                         
            xiRebuild_k = zeros(size(U,2),1);   % reset local xi to 0 as we rebuild the ROM
            
            dr = reduced_constant_vector(d,V,3);
    
            % solve EoMs to get updated nominal solutions eta and dot{eta} (on the deformed mesh
            tic 
            fprintf('____________________\n')
            fprintf('Solving EoMs and sensitivity...\n') 
            TI_NL_PROM = solve_EoMs_and_sensitivities(V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom,h,tmax);                        
            toc
            
                
            eta_0k = TI_NL_PROM.Solution.q;
            eta_k = TI_NL_PROM.Solution.q;
            etad_k = TI_NL_PROM.Solution.qd;
            etadd_k = TI_NL_PROM.Solution.qdd;
            S = TI_NL_PROM.Solution.s;
            Sd = TI_NL_PROM.Solution.sd;
            Sdd = TI_NL_PROM.Solution.sdd;
    
            N = size(eta_k,2);

            uTail = zeros(3,tmax/h);
            for a=1:tmax/h
                uTail(:,a) = V(tailProperties.tailNode*3-2:tailProperties.tailNode*3,:)*TI_NL_PROM.Solution.q(:,a);
            end 
            subplot(2,1,1);
            plot(timePlot,x0Tail+uTail(1,:),'DisplayName',strcat('k=',num2str(k)))
            legend
            drawnow
     
            subplot(2,1,2);
            plot(timePlot,uTail(2,:),'DisplayName',strcat('k=',num2str(k)))
            legend
            drawnow
            
            % solve sensitivity equation 
            % tic
            % fprintf('____________________\n')
            % fprintf('Solving sensitivity...\n') 
            % TI_sens = solve_sensitivities(V,xiRebuild_k,PROM_Assembly, ...
            %    tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom, ...
            %    TI_NL_PROM.Solution.q,TI_NL_PROM.Solution.qd,TI_NL_PROM.Solution.qdd, ...
            %    h,tmax);
            % toc 
            % 
            % S = TI_sens.Solution.q;
            % Sd = TI_sens.Solution.qd;
            % Sdd = TI_sens.Solution.qdd;
        else
            % approximate new solution under new xi, using sensitivity
            fprintf('____________________\n')
            fprintf('Approximating solutions...\n')
            if size(xi_k,1)>1
                S=tensor(S);
                eta_k = eta_0k + double(ttv(S,xiRebuild_k,2));
            else
                eta_k = eta_0k + S*xiRebuild_k;
            end
            
        end 

        % compute cost function and its gradient
        fprintf('____________________\n')
        fprintf('Computing cost function and its gradient...\n') 
        nablaLr = gradient_cost_function_w_constraints_TET4(dr,xiRebuild_k,xi_k,x0,eta_k,etad_k,etadd_k,S,Sd,Sdd,tailProperties,spineProperties,A,b,barrierParam,V);
        LrEvo = [LrEvo, reduced_cost_function_w_constraints_TET4(N,tailProperties,spineProperties,x0,eta_k,etad_k,etadd_k,xiRebuild_k,xi_k,dr,A,b,barrierParam,V)];
        nablaEvo = [nablaEvo,nablaLr];

        % update optimal parameter
        fprintf('____________________\n')
        fprintf('Updating optimal parameter...\n') 
        gradientWeights = adapt_learning_rate(nablaEvo);

        xi_k = xi_k - gStepSize*diag(gradientWeights)*nablaLr
        xiRebuild_k = xiRebuild_k - gStepSize*diag(gradientWeights)*nablaLr;

        xiEvo = [xiEvo,xi_k];
        
        % possible exit conditions
        if size(xi_k,1) >1
            if norm(xiEvo(:,end)-xiEvo(:,end-1))<convCrit
                fprintf('Convergence criterion of %.2g (parameters) fulfilled\n',convCrit)
                break
            elseif length(LrEvo)>5
                if norm(LrEvo(end) - mean(LrEvo(end-4:end))) < convCritCost
                    fprintf('Convergence criterion of %.2g (cost) fulfilled\n',convCrit*100)
                    break
                end
            end
        else
            if norm(xiEvo(end)-xiEvo(end-1))<convCrit
                fprintf('Convergence criterion of %.2g (parameters) fulfilled\n',convCrit)
                break
            elseif length(LrEvo)>5
                if norm(LrEvo(end) - mean(LrEvo(end-4:end))) < convCritCost
                    fprintf('Convergence criterion of %.2g (cost) fulfilled\n',convCrit*100)
                    break
                end
            end
        end

        if k == maxIteration
            fprintf('Maximum number of %d iterations reached\n',maxIteration)
        end
    end
    
    xiStar = xi_k;

end

% Parse input _____________________________________________________________
function [maxIteration,convCrit,convCritCost,barrierParam,gStepSize,nRebuild,rebuildThreshold,FORMULATION,VOLUME,USEJULIA] = parse_inputs(varargin)
    defaultMaxIteration = 50;
    defaultConvCrit = 0.001;
    defaultConvCritCost = 0.1;
    defaultBarrierParam = 500;
    defaultGStepSize = 0.1;
    defaultNRebuild = 10;
    defaultRebuildThreshold = 0.2;
    defaultFORMULATION = 'N1';
    defaultVOLUME = 1;
    defaultUSEJULIA = 0; 
    p = inputParser;
    addParameter(p,'maxIteration',defaultMaxIteration, @(x)validateattributes(x, ...
                    {'numeric'},{'nonempty','integer','positive'}) );
    addParameter(p,'convCrit',defaultConvCrit,@(x)validateattributes(x, ...
                    {'numeric'},{'nonempty','positive'}) );
    addParameter(p,'convCritCost',defaultConvCritCost,@(x)validateattributes(x, ...
                    {'numeric'},{'nonempty','positive'}) );
    addParameter(p,'barrierParam',defaultBarrierParam,@(x)validateattributes(x, ...
                    {'numeric'},{'nonempty','positive'}) );
    addParameter(p,'gStepSize',defaultGStepSize,@(x)validateattributes(x, ...
                    {'numeric'},{'nonempty','positive'}) );
    addParameter(p,'nRebuild',defaultNRebuild,@(x)validateattributes(x, ...
                    {'numeric'},{'nonempty','positive'}) );
    addParameter(p,'rebuildThreshold',defaultRebuildThreshold,@(x)validateattributes(x, ...
                    {'numeric'},{'nonempty','positive'}) );
    addParameter(p,'FORMULATION',defaultFORMULATION,@(x)validateattributes(x, ...
                    {'char'},{'nonempty'}))
    addParameter(p,'VOLUME',defaultVOLUME,@(x)validateattributes(x, ...
                    {'numeric'},{'nonempty'}) );
    addParameter(p,'USEJULIA',defaultUSEJULIA,@(x)validateattributes(x, ...
                    {'numeric'},{'nonempty'}) );
    
    parse(p,varargin{:});
    
    maxIteration = p.Results.maxIteration;
    convCrit = p.Results.convCrit;
    convCritCost = p.Results.convCritCost;
    barrierParam = p.Results.barrierParam;
    gStepSize = p.Results.gStepSize;
    nRebuild = p.Results.nRebuild;
    rebuildThreshold = p.Results.rebuildThreshold;
    FORMULATION = p.Results.FORMULATION;
    VOLUME = p.Results.VOLUME;
    USEJULIA = p.Results.USEJULIA;
end

% Check condition for rebuild _____________________________________________
function cond = check_cond_rebuild(k,lastRebuild,nRebuild, xiRebuild_k, ...
                                    rebuildThreshold,maxIteration)
    cond = 0;

    if mod(k-lastRebuild,nRebuild) == 0 
        cond = 1;
        fprintf('____________________\n')
        fprintf('Rebuilding PROM (max lin. iterations) ...\n')
    elseif any(xiRebuild_k > rebuildThreshold)
        cond = 1;fprintf('____________________\n')
        fprintf('Rebuilding PROM (xi>threshold) ...\n')
    elseif maxIteration-k<0.25*maxIteration ...
            && mod(k-lastRebuild,int16(nRebuild/2)) == 0
        cond = 1;
        fprintf('____________________\n')
        fprintf('Rebuilding PROM (max lin. iteration - close to end) ...\n')
    end
    
end

% Adapt gradient __________________________________________________________
function gradientWeights = adapt_learning_rate(nablaEvo)
    nParam = size(nablaEvo,1);
    gradientWeights = ones(1,nParam);
    for p = 1:nParam
        if sign(nablaEvo(p,end-1)) ~= sign(nablaEvo(p,end)) ...
                && nablaEvo(p,end-1) ~= 0
            fprintf('')
            gradientWeights(p) = 0.5*gradientWeights(p);
            fprintf('Adapting learning rate (change in gradient sign) - xi%d...\n',p)
        end
    end
end


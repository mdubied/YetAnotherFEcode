classdef ImplicitNewmark < handle
    %% Class for Implicit Newmark time Integration 
    % cf. flowchart 7.21 in the following reference:
    % Geradin, M., & Rixen, D. (2015). Mechanical Vibrations : Theory and Application to Structural Dynamics (Third Edition). 
    properties
        % Time integration parameters
        alpha = 0.005
        beta
        h = 0 % time step size
        gamma
        tol = 1e-6      % Relative error tolerance
        
        
        Solution        % Solution data structure
        MaxNRit = 10
        ATS = false     % where adaptive time stepping should be on (true) or not (false)
        hmin = 0        % minimum timestep size (only used when ATS = true)
        NROpt = 3       % Maximum no. of N-R Iterations
        linear = false  % whether system is linear or not
        combinedSensitivity = false    % solve the sensitivities in combination with the EoMs
    end
    
    methods
        function TI = ImplicitNewmark(varargin)
            %% Input parsing
            p = inputParser;
            addParameter(p,'timestep',TI.h, @(x)validateattributes(x, ...
                {'numeric'},{'nonempty'}) );
            addParameter(p,'alpha',TI.alpha, @(x)validateattributes(x, ...
                {'numeric'},{'nonempty'}) );
            addParameter(p,'RelTol',TI.tol, @(x)validateattributes(x, ...
                {'numeric'},{'nonempty','positive'}) );
            addParameter(p,'MaxNRit',TI.MaxNRit, @(x)validateattributes(x, ...
                {'numeric'},{'nonempty','integer','positive'}) );
            addParameter(p,'linear', TI.linear, @(x)validateattributes(x,{'logical'},{'nonempty'}));
            addParameter(p,'hmin', TI.hmin, @(x)validateattributes(x, ...
                {'numeric'},{'nonempty'}) );
            addParameter(p,'ATS', TI.ATS, @(x)validateattributes(x,{'logical'},{'nonempty'}));
            addParameter(p,'combinedSensitivity', TI.combinedSensitivity, @(x)validateattributes(x,{'logical'},{'nonempty'}));
                        
            parse(p,varargin{:});
            
            %% Properties assignment
            TI.alpha = p.Results.alpha;
            TI.beta = (1 + TI.alpha)^2/4;
            TI.gamma = 1/2 + TI.alpha;
            
            TI.h = p.Results.timestep;
            TI.tol = p.Results.RelTol;
            TI.MaxNRit = p.Results.MaxNRit;
            TI.hmin = p.Results.hmin;
            TI.linear = p.Results.linear;
            TI.combinedSensitivity = p.Results.combinedSensitivity;
        end
        function Integrate(obj,x0,xd0,xdd0,tmax,Residual,varargin)            
            % Integrates with Initial condition x0,xd0 from [0 tmax]
            % Residual is a function handle that has the following syntax
            %% Input parsing
            p = inputParser;
            defaultResidualSens = 0;
            defaultActuOnly = false;
            defaultSens0 = 0;
            addParameter(p,'ResidualSens',defaultResidualSens);
            addParameter(p,'actuOnly',defaultActuOnly,@(x)validateattributes(x,{'logical'},{'nonempty'}));
            addParameter(p,'s0',defaultSens0);
            addParameter(p,'sd0',defaultSens0);
            addParameter(p,'sdd0',defaultSens0);
            parse(p,varargin{:});
            ResidualSens = p.Results.ResidualSens;
            actuOnly = p.Results.actuOnly;
            s0 = p.Results.s0;
            sd0 = p.Results.sd0;
            sdd0 = p.Results.sdd0;

            %% Initialize
            if obj.h ==0
                error('Please specify a positive time step')
            end
            
            tic
            t=0;
            time = t;
            q = x0;
            qd = xd0;
            qdd = xdd0;
            q_old = x0;
            qd_old = xd0;
            qdd_old = xdd0;
            NR = 0;
            R = 0;
            i = 1;

            % sensitivity if solved in the combined set up
            if obj.combinedSensitivity
                s = s0;
                sd = sd0;
                sdd = sdd0;
                s_old = s0;
                sd_old = sd0;
                sdd_old = sdd0;
            end
            
            while t < tmax
                t = t+obj.h;
                i = i+1;
                [q_new,qd_new,qdd_new] = obj.Prediction(q_old,qd_old,qdd_old);                
                
                it = -1; % iteration counter
                
                %% linear case
                if obj.linear 
                    it = it + 1; 
                    [r, drdqdd, drdqd, drdq] = Residual(q_new,qd_new,qdd_new,t);
                    S = drdqdd + obj.gamma * obj.h * drdqd + obj.beta * obj.h^2 * drdq;
                    Da = -S\r;
                    [q_new,qd_new,qdd_new] = obj.Correction(q_new,qd_new,qdd_new,Da);
                    epsilon = 0;
                %% Nonlinear case    
                else 
                    %% Newton-Raphson iterations
                    while true                         
                        it = it + 1;
                        
                        %% Compute Residual and Tangent operators
                        [r, drdqdd, drdqd, drdq, c0] = Residual(q_new,qd_new,qdd_new,t);                        
                        
                        %% Check convergence
                        epsilon = norm(r)/c0;
                        if it>cast(0.9*obj.MaxNRit,'int16')
                            disp(['Iteration ' num2str(it) ', Residual norm = '  num2str(epsilon)])
                        end
                        if (epsilon<obj.tol)  % Error < Tolerance : break
                            break;
                        else % Error >= Tolerance : perform correction
                            S = drdqdd + obj.gamma * obj.h * drdqd + obj.beta * obj.h^2 * drdq;
                            Da = -S\r;
                            [q_new,qd_new,qdd_new] = obj.Correction(q_new,qd_new,qdd_new,Da);
                        end
                        
                        %% Adapt time step to maintain an optimal number (obj.NROpt) of N-R iterations 
                        if obj.h > obj.hmin && obj.ATS
                            obj.h = max(obj.h*obj.NROpt/it, obj.hmin);
                        end
                        
                        %% When too many iterations
                        if (it > obj.MaxNRit)
                            warning('Max N-R iterations reached')                            
                            if  epsilon > 1 
                                disp('Exiting time integration: Too high a residual')
                                soltime=toc;
                                obj.Solution.time = time;
                                obj.Solution.q = q;
                                obj.Solution.qd = qd;
                                obj.Solution.qdd = qdd;
                                obj.Solution.NR = NR;
                                obj.Solution.R = R;
                                obj.Solution.soltime = soltime;
                                return
                            else
                                disp('Continuing iterations anyway since the residual is low')                            
                            end
                        end                    

                    end
                end
                
                %% Update solution
                time = [time t];
                NR = [NR it];
                R = [R epsilon];
                if mod(100* t/tmax,20) < 0.35
                    disp(['time integration completed: ', num2str(100* t/tmax), '%'])
                end
%                 disp(['time integration completed: ', num2str(100* t/tmax), '%'])
                
                % needed when solving sensitivity as a separated linear
                % problem
                if size(q_new,2)>1
                    q = cat(3,q,q_new);
                    qd = cat(3,qd,qd_new);
                    qdd = cat(3,qdd,qdd_new);
                    q_old = q_new;
                    qd_old = qd_new;
                    qdd_old = qdd_new;
                else
                    q = [q q_new];
                    qd = [qd qd_new];
                    qdd = [qdd qdd_new];
                    q_old = q_new;
                    qd_old = qd_new;
                    qdd_old = qdd_new;
                end

                % solve the sensitivity as a combined problem
                if obj.combinedSensitivity
                    [s_new,sd_new,sdd_new] = obj.Prediction(s_old,sd_old,sdd_old); 
                    if actuOnly  
                        rSens = ResidualSens(s_new,sd_new,sdd_new,t,q_new,drdqdd, drdqd, drdq);
                    else
                        rSens = ResidualSens(s_new,sd_new,sdd_new,q_new,qd_new,qdd_new,drdqdd, drdqd, drdq,t);
                    end

                    % use the same Jacobian as the one for the EoMs
                    deltaS = -S\rSens;
                    [s_new,sd_new,sdd_new] = obj.Correction(s_new,sd_new,sdd_new,deltaS);

                    if size(s_new,2)>1
                        s = cat(3,s,s_new);
                        sd = cat(3,sd,sd_new);
                        sdd = cat(3,sdd,sdd_new);
                        s_old = s_new;
                        sd_old = sd_new;
                        sdd_old = sdd_new;
                    else
                        s = [s s_new];
                        sd = [sd sd_new];
                        sdd = [sdd sdd_new];
                        s_old = s_new;
                        sd_old = sd_new;
                        sdd_old = sdd_new;
                    end
                    
                end


                
                
            end
            soltime = toc;
            obj.Solution.time = time;
            obj.Solution.q = q;
            obj.Solution.qd = qd;
            obj.Solution.qdd = qdd; 
            obj.Solution.NR = NR;
            obj.Solution.R = R;
            obj.Solution.soltime = soltime;

            if obj.combinedSensitivity
                obj.Solution.s = s;
                obj.Solution.sd = sd;
                obj.Solution.sdd = sdd;
            end
        end
        function[q,qd,qdd] = Prediction(obj,q0,qd0,qdd0)
            qd = qd0 + obj.h * (1 - obj.gamma) * qdd0;
            q = q0 + obj.h * qd0 + (0.5-obj.beta) * obj.h^2 * qdd0;
            qdd = zeros(size(q0));
        end
        function [q,qd,qdd] = Correction(obj,q,qd,qdd,Da)
            q = q + obj.beta * obj.h^2 * Da;
            qd = qd + obj.gamma * obj.h * Da;
            qdd = qdd + Da;
        end
    end
end

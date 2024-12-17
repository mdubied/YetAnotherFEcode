classdef RungeKutta4 < handle
    properties
        h           % Time step size
        Solution    % Structure to store solution data
    end
    
    methods
        function obj = RungeKutta4(varargin)
            % Parse input for time step
            p = inputParser;
            addParameter(p, 'timestep', 0.01, @(x) validateattributes(x, {'numeric'}, {'positive'}));
            parse(p, varargin{:});
            
            % Assign time step
            obj.h = p.Results.timestep;
        end
        
        function Integrate(obj, y0, tmax, Dynamics)
            % Initialize solution data
            t = 0;
            time = t;
            y = y0;
            
            % Solution data structure
            obj.Solution.time = time;
            obj.Solution.y = y;
            
            % Runge-Kutta 4 integration loop
            while t < tmax
                t = t + obj.h;
                
                % Compute RK4 terms
                k1 = Dynamics(y, t);
                k2 = Dynamics(y + 0.5 * obj.h * k1, t + 0.5 * obj.h);
                k3 = Dynamics(y + 0.5 * obj.h * k2, t + 0.5 * obj.h);
                k4 = Dynamics(y + obj.h * k3, t + obj.h);
                
                % Update solution
                y = y + (obj.h / 6) * (k1 + 2*k2 + 2*k3 + k4);
                
                % Store results
                time = [time; t];
                obj.Solution.y = [obj.Solution.y, y];
                obj.Solution.time = time;
                
                if mod(100* t/tmax,2) < 0.35
                    disp(['time integration completed: ', num2str(100* t/tmax), '%'])
                end
            end
        end
    end
end

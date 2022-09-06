classdef Tri3Element < ContinuumElement
    properties
        nodes = []          % global coordinates of element nodes
        nodeIDs = []        % the index location of element nodes
        nDOFPerNode = 2     % number of DOFs per node
        nNodes = 3          % number of nodes per element
        nDim = 2            % number of dimensions in local coordinates
        nelDOFs
        quadrature       	% weights and points for gauss quadrature
        Material           	% Object of class Material
        initialization     	% some 0-matrices to speedup numerical integration
        elType = 'TRI'
    end
    
    properties
        thickness = 0       % element thickness, by default zero    
    end
    
    methods
        
        % MINIMUM REQUIRED FUNCTIONS ______________________________________
        
        function self = Tri3Element(thickness, Material, Ngauss)
            % Self function (constructor)
            if nargin == 2
                Ngauss = 2;
            end
            self.thickness = thickness;
            self.nelDOFs = self.nNodes * self.nDOFPerNode;
            ContinuumElementConstructor(self, Material, Ngauss);
        end
        
        function [G,detJ,dH] = shape_function_derivatives(self, X)
            %______________________________________________________________
            %
            % [G,detJ,dH] = shape_function_derivatives(self, X)
            % G = shape function derivatives in physical coordinates, s.t.
            % th=G*p with th={ux uy vx vy}' (ux=du/dx...)
            % and p={u1,v1,...,u4,v4}' (nodal displacements).
            % detJ = det(J), J=jacobian
            %______________________________________________________________
            r = X(1);
            s = X(2);
            xy = self.nodes;
            % shape function derivatives in natural coordinates
            dHn = [ -1,  -1; 
                     1,  0; 
                     0,  1].';
            J = dHn*xy;
            detJ = det(J);
            dH = J\dHn;	% derivatives in physical coordinates,
                      	% 2x8 matrix, [dNi_dx; dNi_dy]
                       	% with i = 1...8
            G = self.initialization.G;
            G(1:2,1:2:end) = dH;
            G(3:4,2:2:end) = dH;
        end
        
        % ADDITIONAL FUNCTIONS ____________________________________________

        function F = force_length_prop_skin_normal(self,specificFace)
            % _____________________________________________________________
            % Apply a force proportional to the the length of a skin
            % element on the 2 corresponding nodes, in the direction normal
            % to the face considered. Robust to elements having up to 2 
            % skin faces.
            % _____________________________________________________________
            F = sparse(self.nelDOFs,1);
            [startNode, endNode, nextNode] = getNodesFromFace(specificFace(1));
            n = normalVector(self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));

            length = norm(self.nodes(endNode,:)-self.nodes(startNode,:));
            F(startNode*2-1) = length/2 * n(1);  % x-coordinate
            F(startNode*2) = length/2 * n(2);    % y-coordinate
            F(endNode*2-1) = length/2 * n(1);    % x-coordinate
            F(endNode*2) = length/2 * n(2);      % y-coordinate
        
            if specificFace(2) ~= 0
                [startNode, endNode, nextNode] = getNodesFromFace(specificFace(2));
                n = normalVector(self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));

                length = norm(self.nodes(endNode,:)-self.nodes(startNode,:));
                F(startNode*2-1) = F(startNode*2-1) + length/2 * n(1);    % x-coordinate
                F(startNode*2) = F(startNode*2) + length/2 * n(2);        % y-coordinate
                F(endNode*2-1) = F(endNode*2-1) + length/2 * n(1);        % x-coordinate
                F(endNode*2) = F(endNode*2) + length/2 * n(2);            % y-coordinate
            end

            function n = normalVector(startNodePos, endNodePos, nextNodePos)
                % Returns a vector of unit length normal to the surface
                % defined by `startNodePos' and enNodePos', pointing
                % outward the considered element
                n = [0 1; -1 0]*[endNodePos(1)-startNodePos(1); endNodePos(2)-startNodePos(2)];  % 90° rotation clockwise
                n = n/norm(n); % normalizing n
                dotProd = dot(n,[nextNodePos(1)-endNodePos(1); nextNodePos(2)-endNodePos(2)]);
                if dotProd >= 0
                    n = -n; % inverting the direction of n if it points toward the element
                end
            end

            function [startNode, endNode, nextNode] = getNodesFromFace(face)
                % Returns the nodes' indexes from a give face. The
                % ordering is done counterclokewise. `nextNode is the index
                % of the node on the next surface.
                if face ~= 3
                    startNode = face;
                    endNode = face + 1;
                else
                    startNode = face;
                    endNode = 1;
                end
                if face ~= 2
                    nextNode = endNode +1;
                else
                    nextNode = 1;
                end
            end 
        end 

        function F = drag_force(self, specificFace, vwater, rho)
            % _____________________________________________________________
            % Returns the drag force F_drag acting the face
            % ``specificFace'' of an element. The force is divided by 2 and
            % applied to the nodes at of the considered face.
            % _____________________________________________________________
            F = sparse(self.nelDOFs,1);
            [startNode, endNode, nextNode] = getNodesFromFace(specificFace(1));

            % Terms needed for the drag force - REPLACE VREL (to do) !!
            n = normalVector(self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));  
            length = norm(self.nodes(endNode,:)-self.nodes(startNode,:));
            A = length; % in 2D, area proportional force is length proportional
            ns = 2; % 2D surfaces of TRI3 are composed of 2 elements
            vrel = vwater - 1/ns*[0;0]; %(self.nodes(startNode,:) + self.nodes(endNode,:)); % need to be replaced by velocity
            d = vrel/norm(vrel);
            Cd = 2*(n.'*d)^2;
            
            % Drag force
            Fdrag = 1/2*rho*Cd*norm(vrel)^2*d;

            % Apply drag force on nodes
            F(startNode*2-1) = F(startNode*2-1) + Fdrag(1)/2;    % x-coordinate
            F(startNode*2) = F(startNode*2) + Fdrag(2)/2;        % y-coordinate
            F(endNode*2-1) = F(endNode*2-1) + Fdrag(1)/2;        % x-coordinate
            F(endNode*2) = F(endNode*2) + Fdrag(2)/2;            % y-coordinate

            if specificFace(2) ~= 0
                [startNode, endNode, nextNode] = getNodesFromFace(specificFace(2));

                % Terms needed for the drag force - REPLACE VREL (to do) !!
                n = normalVector(self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));  
                length = norm(self.nodes(endNode,:)-self.nodes(startNode,:));
                A = length; % in 2D, area proportional force is length proportional
                ns = 2; % 2D surfaces of TRI3 are composed of 2 elements
                vrel = vwater - 1/ns*[0;0]; %(self.nodes(startNode,:) + self.nodes(endNode,:)); % need to be replaced by velocity
                d = vrel/norm(vrel);
                Cd = 2*(n.'*d)^2;
                
                % Drag force
                Fdrag = 1/2*rho*Cd*norm(vrel)^2*d;
    
                % Apply drag force on nodes
                F(startNode*2-1) = F(startNode*2-1) + Fdrag(1)/2;    % x-coordinate
                F(startNode*2) = F(startNode*2) + Fdrag(2)/2;        % y-coordinate
                F(endNode*2-1) = F(endNode*2-1) + Fdrag(1)/2;        % x-coordinate
                F(endNode*2) = F(endNode*2) + Fdrag(2)/2;            % y-coordinate
            end
            
            function n = normalVector(startNodePos, endNodePos, nextNodePos)
                % Returns a vector of unit length normal to the surface
                % defined by `startNodePos' and enNodePos', pointing
                % outward the considered element
                n = [0 1; -1 0]*[endNodePos(1)-startNodePos(1); endNodePos(2)-startNodePos(2)];  % 90° rotation clockwise
                n = n/norm(n); % normalizing n
                dotProd = dot(n,[nextNodePos(1)-endNodePos(1); nextNodePos(2)-endNodePos(2)]);
                if dotProd >= 0
                    n = -n; % inverting the direction of n if it points toward the element
                end
            end

            function [startNode, endNode, nextNode] = getNodesFromFace(face)
                % Returns the nodes' indexes from a give face. The
                % ordering is done counterclokewise. `nextNode is the index
                % of the node on the next surface.
                if face ~= 3
                    startNode = face;
                    endNode = face + 1;
                else
                    startNode = face;
                    endNode = 1;
                end
                if face ~= 2
                    nextNode = endNode +1;
                else
                    nextNode = 1;
                end
            end 

        
        end


    end 
        
    
    methods (Static)
        
        function N = shape_functions(X)
            % N = shape_functions(g,h,r)
            % SHAPE FUNCTIONS FOR A 3-NODED TRIANGLE
            % - see Abaqus documentation:
            % Abaqus theory guide > Elements > Continuum elements > ...
            % ... > 3.2.6 Triangular, tetrahedral, and wedge elements
            g = X(1);
            h = X(2);
            N = [ (1 - g - h);g ;h ];                            
        end
        
        function X = natural_coordinates
            X = [0  0       % node 1
                 1  0       % node 2
                 0  1];     % node 3
        end
        
    end
        
end % classdef
    

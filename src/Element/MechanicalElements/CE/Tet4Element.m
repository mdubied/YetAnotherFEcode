classdef Tet4Element < ContinuumElement
    properties
        nodes = []          % global coordinates of element nodes
        nodeIDs = []        % the index location of element nodes
        nDOFPerNode = 3     % number of DOFs per node
        nNodes = 4          % number of nodes per element
        nDim = 3            % number of dimensions in local coordinates
        nelDOFs
        quadrature       	% weights and points for gauss quadrature
        Material           	% Object of class Material
        initialization     	% some 0-matrices to speedup numerical integration
        elType = 'TET'
    end
    
    properties
        thickness = 1       % element thickness  
    end
    
    methods
        
        % MINIMUM REQUIRED FUNCTIONS ______________________________________
        
        function self = Tet4Element(Material)
            % _____________________________________________________________
            %
            % SELF-FUNCTION
            % self = Tet4Element(Material,Ngauss)
            % defines element's properties
            %______________________________________________________________
            Ngauss = 1; % note: shape function derivatives are constants,
                        % and M is a lumped-parameter mass matrix. No
                        % quadrature integration is actually needed.
            self.thickness = 1;
            self.nelDOFs = self.nNodes * self.nDOFPerNode;
            ContinuumElementConstructor(self, Material, Ngauss);
        end
        
        function Mel = mass_matrix(self)
            % _____________________________________________________________
            %
            % Mel = mass_matrix_global(self,nodes,~);
            % Mel: element-level mass matrix (in global coordinates)
            %______________________________________________________________
            rho = self.Material.DENSITY;
            m = self.vol * rho / 4; % lumped masses are better for TET4
            Mel = sparse(eye(12)*m);
        end
        
        function [G,detJ,dH] = shape_function_derivatives(self, X)
            %______________________________________________________________
            %
            % [G,detJ,dH] = shape_function_derivatives(self)
            % G = shape function derivatives in physical coordinates, s.t.
            % th=G*p with th={ux uy uz vx vy vz wx wy wz}' (ux=du/dx...)
            % and p={u1,v1,w1,...,u4,v4,w4}' (nodal displacements).
            % detJ = det(J), J=jacobian
            %______________________________________________________________
            xyz = self.nodes;
            % First order thetrahedron. Shape functions:
            %   N = [(1-g-h-r), g, h, r].';
            % Shape function derivatives in natural coordinates:
            dHn = [ -1, 1, 0, 0;
                    -1, 0, 1, 0;
                    -1, 0, 0, 1];    
            
            J = dHn*xyz;
            J1 = [0 0 0; 0 0 0; 0 0 0];
            J1(1,1) = J(2,2)*J(3,3) - J(2,3)*J(3,2);
            J1(2,1) = J(2,3)*J(3,1) - J(2,1)*J(3,3);
            J1(3,1) = J(2,1)*J(3,2) - J(2,2)*J(3,1);
            J1(1,2) = J(1,3)*J(3,2) - J(1,2)*J(3,3);
            J1(2,2) = J(1,1)*J(3,3) - J(1,3)*J(3,1);
            J1(3,2) = J(1,2)*J(3,1) - J(1,1)*J(3,2);
            J1(1,3) = J(1,2)*J(2,3) - J(1,3)*J(2,2);
            J1(2,3) = J(1,3)*J(2,1) - J(1,1)*J(2,3);
            J1(3,3) = J(1,1)*J(2,2) - J(1,2)*J(2,1);
            detJ = J(1,1)*J1(1,1) + J(1,2)*J1(2,1) + J(1,3)*J1(3,1);
            J1 = J1/detJ;
            dH = J1*dHn;   	% derivatives in physical coordinates,
                            % 3x4 matrix, [dNi_dx; dNi_dy; dNi_dz]
                            % with i = 1...4
            G = self.initialization.G;
            G(1:3,1:3:12) = dH;
            G(4:6,2:3:12) = dH;
            G(7:9,3:3:12) = dH;
        end
        
        % HYDRODYNAMIC FORCES _____________________________________________

        function T = spine_momentum_tensor(self, spineNodeIndexInElement, normalisation)
            % _____________________________________________________________ 
            % Compute the tensor related to the spine change in momentum
            % _____________________________________________________________
           
            % get matrix corresponding to the configuration
            A = A_TET4(spineNodeIndexInElement(1));
            B = B_TET4(spineNodeIndexInElement(1),spineNodeIndexInElement(2));
            R = [0 -1 0 0 0 0 0 0 0 0 0 0;
                 1 0 0 0 0 0 0 0 0 0 0 0;
                 0 0 0 0 0 0 0 0 0 0 0 0;
                 0 0 0 0 -1 0 0 0 0 0 0 0;
                 0 0 0 1 0 0 0 0 0 0 0 0;
                 0 0 0 0 0 0 0 0 0 0 0 0;
                 0 0 0 0 0 0 0 -1 0 0 0 0;
                 0 0 0 0 0 0 1 0 0 0 0 0;
                 0 0 0 0 0 0 0 0 0 0 0 0;
                 0 0 0 0 0 0 0 0 0 0 -1 0;
                 0 0 0 0 0 0 0 0 0 1 0 0;
                 0 0 0 0 0 0 0 0 0 0 0 0];     % 90 degrees rotation counterclock-wise around z axis
            
            % compute force
            T = -0.25*pi*normalisation*einsum('mJ,mi,iK,Is,sL->IJKL',A,R,B,R,B);
            T = double (T);

        end

        % Test functions
        function F = force_length_prop_skin_normal(self,specificFace)
            % _____________________________________________________________
            % Applies a force proportional to the the area of a skin
            % element on the 3 corresponding nodes, in the direction normal
            % to the face considered. Robust to elements having up to 2 
            % skin faces.
            % _____________________________________________________________
            F = sparse(self.nelDOFs,1);
            [startNode, midNode, endNode, nextNode] = get_node_from_face(self, specificFace(1));
            n = normal_vector(self, self.nodes(startNode,:), self.nodes(midNode,:),self.nodes(endNode,:), self.nodes(nextNode,:));

            area = norm(cross(self.nodes(midNode,:)-self.nodes(startNode,:),self.nodes(endNode,:)-self.nodes(startNode,:)));
            Force = area*n;
            F = apply_force(self, F, Force, startNode, midNode, endNode);
        
            if specificFace(2) ~= 0
                [startNode, midNode, endNode, nextNode] = get_node_from_face(self, specificFace(2));
                n = normal_vector(self, self.nodes(startNode,:), self.nodes(midNode,:),self.nodes(endNode,:), self.nodes(nextNode,:));
    
                area = norm(cross(self.nodes(midNode,:)-self.nodes(startNode,:),self.nodes(endNode,:)-self.nodes(startNode,:)));
                Force = area*n;
                F = apply_force(self, F, Force, startNode, midNode, endNode);
            end
        end 

        function F = apply_force(~, Fbase, FtoApply, startNode, midNode, endNode)
            % _____________________________________________________________
            % Returns the addition of the force to apply FtoApply and the
            % base force Fbase.
            % _____________________________________________________________
            F = Fbase;
            F(startNode*3-2) = Fbase(startNode*3-2) + FtoApply(1)/3;    % x-coordinate
            F(startNode*3-1) = Fbase(startNode*3-1) + FtoApply(2)/3;    % y-coordinate
            F(startNode*3)   = Fbase(startNode*3) + FtoApply(3)/3;      % z-coordinate
            F(midNode*3-2)   = Fbase(midNode*3-2) + FtoApply(1)/3;      % x-coordinate
            F(midNode*3-1)   = Fbase(midNode*3-1) + FtoApply(2)/3;      % y-coordinate
            F(midNode*3)     = Fbase(midNode*3) + FtoApply(3)/3;        % z-coordinate
            F(endNode*3-2)   = Fbase(endNode*3-2) + FtoApply(1)/3;      % x-coordinate
            F(endNode*3-1)   = Fbase(endNode*3-1) + FtoApply(2)/3;      % y-coordinate
            F(endNode*3)     = Fbase(endNode*3) + FtoApply(3)/3;        % z-coordinate
        end 
        
        % Useful helper functions
        function n = normal_vector(~, startNodePos, midNodePos, endNodePos, nextNodePos)
            % _____________________________________________________________
            % Returns a vector of unit length normal to the surface
            % defined by `startNodePos', 'midNodePos', and 'enNodePos',
            % pointing outward the considered element. 'nextNodePos' is the
            % position of last node of the element.
            % _____________________________________________________________
            v1 = midNodePos - startNodePos;
            v2 = endNodePos - startNodePos;
            v3 = nextNodePos - startNodePos;
            n = cross(v1,v2);
            n = n/norm(n); % normalizing n
            dotProd = dot(n,v3);
            if dotProd >= 0
                n = -n; % inverting the direction of n if it points toward the element
                disp("inversion")
            else
                disp("no inversion")
            end
        end


        function inv = normal_vector_inversion(~, startNodePos, midNodePos, endNodePos, nextNodePos)
            % _____________________________________________________________
            % Returns returns -1 if the normal vector computed with the
            % basic formula cross(v1,v2)/norm(cross(v1,v2)) (see
            % `normal_vector') needs to be inverted. Returns 1 if this is
            % not the case.
            % _____________________________________________________________
            v1 = midNodePos - startNodePos;
            v2 = endNodePos - startNodePos;
            v3 = nextNodePos - startNodePos;
            n = cross(v1,v2);
            n = n/norm(n); % normalizing n
            dotProd = dot(n,v3);
            if dotProd >= 0
                inv = -1; % inversion is needed if n points toward the element
            else
                inv = 1;
            end
        end


        function [startNode, midNode, endNode, nextNode] = get_node_from_face(~, face)
            % _____________________________________________________________
            % Returns the nodes' indexes from a give face. 
            % `nextNode' is the index off the 4th element's node that is 
            % not part of the face.
            % _____________________________________________________________
            switch face
                case 1
                    startNode = 1;
                    midNode = 2;
                    endNode = 3;
                    nextNode = 4;
                case 2
                    startNode = 2;
                    midNode = 3;
                    endNode = 4;
                    nextNode = 1;
                case 3
                    startNode = 3;
                    midNode = 4;
                    endNode = 1;
                    nextNode = 2;
                case 4
                    startNode = 4;
                    midNode = 1;
                    endNode = 2;
                    nextNode = 3;
                otherwise
                    disp("Error in `get_node_from_face', Tet4Element.m")
            end

            
        end
    end
    
    methods (Static)
        
        function N = shape_functions(X)
            % N = shape_functions(g,h,r)
            % SHAPE FUNCTIONS FOR A 4-NODED TETRAHEDRON
            % - see Abaqus documentation:
            % Abaqus theory guide > Elements > Continuum elements > ...
            % ... > 3.2.6 Triangular, tetrahedral, and wedge elements
            g = X(1);
            h = X(2);
            r = X(3);
            N = [(2*(1-g-h-r)-1)*(1-g-h-r)
                    (2*g-1)*g
                    (2*h-1)*h
                    (2*r-1)*r];
        end
        
        function X = natural_coordinates
            X = [ ...
                0 0 0
                1 0 0
                0 1 0
                0 0 1];
        end
        
    end

        
end % classdef
    

classdef Quad8Element < Element
    properties
        nodes = []          % global coordinates of element nodes
        nodeIDs = []        % the index location of element nodes
        nDOFPerNode = 2     % number of DOFs per node
        nNodes = 8          % number of nodes per element
        nDim = 2            % number of dimensions in local coordinates
    end
    
    properties
        thickness = 0	% element thickness, by default zero
        Material       	% Object of class Material
        quadrature    	% weights and points for gauss quadrature
        initialization 	% some 0-matrices to speedup numerical integration     
    end
    
    properties (Dependent)
        uniformBodyForce
        area                % area of the element
    end
    
    methods
        
        % MINIMUM REQUIRED FUNCTIONS ______________________________________
        
        function self = Quad8Element(thickness, Material, Ngauss)
            % _____________________________________________________________
            %
            % SELF-FUNCTION
            % self = Quad8Element(Material,Ngauss)
            % defines element's properties
            %______________________________________________________________
            self.Material = Material;
            if nargin == 2
                Ngauss = 2;
            end
            [x,w]=lgwt(Ngauss,-1,1);
            X = zeros(2,Ngauss^2);
            W = zeros(Ngauss^2,1);
            cont = 1;
            for ii = 1:Ngauss
                for jj = 1:Ngauss
                    X(:,cont) = [x(ii) x(jj)].';
                    W(cont) = w(ii)*w(jj);
                    cont = cont+1;
                end
            end
            self.quadrature.Ng = Ngauss;
            self.quadrature.X = X;      % gauss integration points
            self.quadrature.W = W;      % gauss integration weights
            self.thickness = thickness;
            % INIZIALIZATION of some matrices (this should speedup
            % numerical integration)
            C = self.Material.get_stress_strain_matrix_2D * thickness;
            H = [1 0 0 0; 
                0 0 0 1; 
                0 1 1 0];
            self.initialization.A = zeros(3,4); % nonlinear strain
            self.initialization.G = zeros(4,16);% shape function derivatives
            self.initialization.Z = zeros(8);   % zero-matrix
            self.initialization.K = zeros(16);  % stiffness-element matrix
            self.initialization.F = zeros(16,1);% internal forces (element)
            self.initialization.C = C;          % constitutive law matrix
            self.initialization.H = H;          % linear strain
        end
        
        function Mel = mass_matrix(self)
            % _____________________________________________________________
            %
            % Mel = mass_matrix_global(self,nodes,~);
            % Mel: element-level mass matrix (in global coordinates)
            %______________________________________________________________
            X = self.quadrature.X;
            W = self.quadrature.W;
            rho = self.Material.DENSITY;
            t = self.thickness;
            Mel = zeros(16);
            for ii = 1:length(W)
                g = X(1,ii);
                h = X(2,ii);
                we = W(ii); % weights
                % shape functions and detJ (ABAQUS ORDER)
                N = self.shape_functions(g,h);
                [~,detJ] = shape_function_derivatives(self,g,h);
                NN = kron(N', eye(2));
                % integration of K and M through GAUSS QUADRATURE
                Mel = Mel + ( NN' * NN )*( we * detJ );
            end
            Mel = sparse(t*rho*Mel);
        end
        
        function [K,F] = tangent_stiffness_and_force(self, x)
            displ = self.extract_element_data(x);
            X = self.quadrature.X;
            W = self.quadrature.W;
            K = self.initialization.K;
            F = self.initialization.F;
            C = self.initialization.C;
            H = self.initialization.H;
            ZZ = self.initialization.Z;
            for ii = 1:length(W)
                g = X(1,ii);
                h = X(2,ii);
                we = W(ii); % weights
                [G,detJ,dH] = shape_function_derivatives(self,g,h);
                th  = G*displ;
                A =	[th(1)	0     th(3) 0;
                     0    	th(2) 0     th(4);
                     th(2)	th(1) th(4) th(3)];
                % Green Strain tensor
                E = (H + 1/2*A)*th;
                % second Piola-Kirchhoff stress tensor
                s = C*E; % s = [S11 S22 S12]
                S = [s(1), s(3); s(3), s(2)];
                Bnl = (H + A)*G;
                % functions to integrate over volume
                int_K1 = Bnl'*C*Bnl;
                HSH = dH'*S*dH;
                int_Ks = [HSH ZZ; ZZ HSH]; % (faster than blkdiag)
                int_K = int_K1 + int_Ks;
                int_F = Bnl'*s;
                % integration of K and F through Gauss quadrature
                K = K + int_K * (we * detJ);
                F = F + int_F * (we * detJ);
            end
        end
        
        function xe = extract_element_data(self, x)
            % x is a vector of full DOFs
            index = get_index(self.nodeIDs, self.nDOFPerNode);
            xe = x(index,:);
        end
        
        function F = uniform_body_force(self,~)
            % _____________________________________________________________
            %
            % F = uniform_body_force(self,direction)
            % This function computes a load along direction=2(Y) by
            % dividing the load on the 8 nodes according to the element
            % area (A/8) [it might not be the best way, but still...]
            %______________________________________________________________
            F = sparse(16,1);
            F(2:2:end) = self.area/8; % uniformly distributed pressure on the structure
        end
        
        function [T2, globalSubs] = T2(self, varargin)
            % this function computes the 3-tensor corresponding to the 
            % quadratic component of the nonlinear internal force in 
            % global coordinates at the element level.
            
            if ~isempty(varargin)
                Ve = varargin{1};
                m = size(Ve, 2);
                Vflag = true;
            else
                m = self.nNodes*self.nDOFPerNode;
                Vflag = false;
            end
            
            % global DOFs associated to the element nodes
            index = get_index(self.nodeIDs,self.nDOFPerNode);
            
            % location of each dimension of tensor in global DOFs
            globalSubs = {index, index, index};
                        
            X = self.quadrature.X;
            W = self.quadrature.W;

            C = self.initialization.C;  % constitutive law matrix
            H = self.initialization.H;  % Linear strain matrix: eps_l = H*th

            % Quadratic strain matrix: A = L.th, eps_quad = A*th
            L = tenzeros([3,4,4]);
            L(1,1,1) = 1; L(3,2,1) = 1; L(3,1,2) = 1; L(2,2,2) = 1;
            L(1,3,3) = 1; L(3,4,3) = 1; L(3,3,4) = 1; L(2,4,4) = 1;
            
            Q3h = tenzeros([m,m,m]);
            for ii = 1:length(W)
                g = X(1,ii);
                h = X(2,ii);
                we = W(ii); % weights
                [G,detJ] = shape_function_derivatives(self,g,h);
                if Vflag
                    G = G*Ve;
                end
                % G(x,y,z) and detJ from the position of the gauss points

                %construct core part of the tensors for each gauss point
                GHC = tensor((C*H*G)');
                TG = tensor(G);  %create tensor object out of matrix
                LGG = ttt(ttt(L,TG,3,1),TG,2,1);

                Q3h_int = ttt(GHC,LGG,2,1);
                Q3h = Q3h + Q3h_int*detJ*we;
            end

            % build third order tensors using Q3h
            Q3ht = permute(Q3h,[3 2 1]);
            T2 = Q3h./2 + Q3ht;           
        end
        
        function [T3, globalSubs] = T3(self, varargin)
            % this function computes the 4-tensor corresponding to the 
            % quadratic component of the nonlinear internal force in 
            % global coordinates at the element level.
            
            if ~isempty(varargin)
                Ve = varargin{1};
                m = size(Ve, 2);
                Vflag = true;
            else
                m = self.nNodes*self.nDOFPerNode;
                Vflag = false;
            end
            
            % global DOFs associated to the element nodes
            index = get_index(self.nodeIDs,self.nDOFPerNode);
            
            % location of each dimension of tensor in global DOFs
            globalSubs = cell(4,1);
            globalSubs(:) = {index};
                        
            X = self.quadrature.X;
            W = self.quadrature.W;

            C = self.initialization.C;  % constitutive law matrix

            % Quadratic strain matrix: A = L.th, eps_quad = A*th
            L = tenzeros([3,4,4]);
            L(1,1,1) = 1; L(3,2,1) = 1; L(3,1,2) = 1; L(2,2,2) = 1;
            L(1,3,3) = 1; L(3,4,3) = 1; L(3,3,4) = 1; L(2,4,4) = 1;
            
            T3 = tenzeros([m,m,m,m]);

            for ii = 1:length(W)
                g = X(1,ii);
                h = X(2,ii);
                we = W(ii); % weights
                [G,detJ] = shape_function_derivatives(self,g,h);
                if Vflag
                    G = G*Ve;
                end
                % G(x,y,z) and detJ from the position of the gauss points

                %construct core part of the tensors for each gauss point
                TC = tensor(C);  %create tensor object, rename it to distinguish
                TG = tensor(G);  %create tensor object out of matrix
                LGG = ttt(ttt(L,TG,3,1),TG,2,1);

                Q4h_int = ttt(ttt(permute(LGG,[2 1 3]),TC,2,1),LGG,3,1);
                T3 = T3 + Q4h_int*detJ*we/2;
            end
           
        end
        
        function [Kd] = stiffness_derivative(self,x)
            % this function computes the element stiffness derivative matrix
            % with respect to the amplitude q of an imposed displacement
            % field U (directional derivative), that is:
            %      dK|         dK(Uq)|
            % Kd = --|       = ----- |
            %      dq|_{q=0}    dq   |_{q=0}
            % and evaluated for q=0.
            % x : element DOF values of U
            
            displ = self.extract_element_data(x);
            X = self.quadrature.X;
            W = self.quadrature.W;
            Kd= self.initialization.K;
            C = self.initialization.C;
            H = self.initialization.H;
            % Quadratic strain matrix: A = L.th, eps_quad = A*th
            L = zeros([3,4,4]);
            L(1,1,1) = 1; L(3,2,1) = 1; L(3,1,2) = 1; L(2,2,2) = 1;
            L(1,3,3) = 1; L(3,4,3) = 1; L(3,3,4) = 1; L(2,4,4) = 1;
            for ii = 1:length(W)
                g = X(1,ii);
                h = X(2,ii);
                we = W(ii); % weights
                [G,detJ] = shape_function_derivatives(self,g,h);
                th  = G*displ;
                A =	[th(1)	0     th(3) 0;
                     0    	th(2) 0     th(4);
                     th(2)	th(1) th(4) th(3)];
                b1 = einsum('ijk,kl->ijl',L,G);
                b2 = einsum('ijk,il->jkl',b1,C*H*th);
                b = G'*b2;
                int_dK = G'*(H'*C*A + A'*C*H)*G + b;
                Kd = Kd + int_dK * detJ * we;
            end
        end
        
        % ANCILLARY FUNCTIONS _____________________________________________
        
        function A = get.area(self)
            % Integrate detJ (jacobian from isoparametric to physical
            % coordinates) over the area to get A
            detJ = 0;
            W = self.quadrature.W;
            X = self.quadrature.X;
            for ii = 1 : length( W )
                for jj = 1 : length( W )
                    g = X(ii);
                    h = X(ii);
                    [~, detJ_i] = shape_function_derivatives(self, g, h);
                    detJ = detJ + detJ_i * W(ii) * W(jj);
                end
            end
            A = detJ;
        end
        
        function [G,detJ,dH] = shape_function_derivatives(self,g,h)
            %______________________________________________________________
            %
            % [G,detJ,dH] = shape_function_derivatives(self,g,h,r)
            % G = shape function derivatives in physical coordinates, s.t.
            % th=G*p with th={ux uy vx vy}' (ux=du/dx...)
            % and p={u1,v1,...,u8,v8}' (nodal displacements).
            % detJ = det(J), J=jacobian
            %______________________________________________________________
            xy = self.nodes;
            % shape function derivatives in natural coordinates
            dHn = [ ...
                    -((h - 1)*(g + h + 1))/4 - ((g - 1)*(h - 1))/4, -((g - 1)*(g + h + 1))/4 - ((g - 1)*(h - 1))/4;
                    ((h - 1)*(h - g + 1))/4 - ((g + 1)*(h - 1))/4,   ((g + 1)*(h - g + 1))/4 + ((g + 1)*(h - 1))/4;
                    ((h + 1)*(g + h - 1))/4 + ((g + 1)*(h + 1))/4,   ((g + 1)*(g + h - 1))/4 + ((g + 1)*(h + 1))/4;
                    ((h + 1)*(g - h + 1))/4 + ((g - 1)*(h + 1))/4,   ((g - 1)*(g - h + 1))/4 - ((g - 1)*(h + 1))/4;
                                                        g*(h - 1),                                     g^2/2 - 1/2;
                                                      1/2 - h^2/2,                                      -h*(g + 1);
                                                       -g*(h + 1),                                     1/2 - g^2/2;
                                                      h^2/2 - 1/2,                                       h*(g - 1)]';
            J = dHn*xy;
            detJ = det(J);
            dH = J\dHn;	% derivatives in physical coordinates,
                      	% 2x16 matrix, [dNi_dx; dNi_dy]
                       	% with i = 1...10
            G = self.initialization.G;
            G(1:2,1:2:end) = dH;
            G(3:4,2:2:end) = dH;
        end
        
    end % methods
    
    methods (Static)
        
        function N = shape_functions(g,h)
            % N = shape_functions(g,h,r)
            % SHAPE FUNCTIONS FOR A 8-NODED QUADRILATERAL
            % - see Abaqus documentation:
            % Abaqus theory guide > Elements > Continuum elements > ...
            % ... > 3.2.4 Solid isoparametric quadrilaterals and hexahedra
            N = 1/4*[...
                    -(1-g)*(1-h)*(1+g+h); 
                    -(1+g)*(1-h)*(1-g+h);
                    -(1+g)*(1+h)*(1-g-h); 
                    -(1-g)*(1+h)*(1+g-h);
                    2*(1-g)*(1+g)*(1-h);  
                    2*(1-h)*(1+h)*(1+g);
                    2*(1-g)*(1+g)*(1+h);  
                    2*(1-h)*(1+h)*(1-g)];
        end
        
        function X = natural_coordinates
            X = [ ...
                -1  -1  % node 1 (corner)
                1   -1  % node 2 (corner)
                1   1	% node 3 (corner)
                -1  1   % node 4 (corner)
                0   -1  % node 5
                1   0   % node 6
                0   1   % node 7
                -1  0]; % node 8
        end
        
    end
        
end % classdef
    

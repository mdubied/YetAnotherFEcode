classdef Hex20Element < Element
    properties
        nodes = []          % global coordinates of element nodes
        nodeIDs = []        % the index location of element nodes
        nDOFPerNode = 3     % number of DOFs per node
        nNodes = 20         % number of nodes per element
        nDim = 3            % number of dimensions in local coordinates
    end
    
    properties
        quadrature              % weights and points for gauss quadrature
        Material                % Object of class Material
        initialization          % some 0-matrices to speedup numerical integration
    end
    
    properties (Dependent)
        uniformBodyForce
        vol                     % volume of the element
    end
    
    methods
        
        % MINIMUM REQUIRED FUNCTIONS ______________________________________
        
        function self = Hex20Element(Material,Ngauss)
            % _____________________________________________________________
            %
            % SELF-FUNCTION
            % self = Hex20Element(Material,Ngauss)
            % defines element's properties
            %______________________________________________________________
            self.Material = Material;
            if nargin == 1
                Ngauss = 2;
            end
            [x,w]=lgwt(Ngauss,-1,1);
            X = zeros(3,Ngauss^3);
            W = zeros(Ngauss^3,1);
            cont = 1;
            for ii = 1:Ngauss
                for jj = 1:Ngauss
                    for kk = 1:Ngauss
                        X(:,cont) = [x(ii) x(jj) x(kk)].';
                        W(cont) = w(ii)*w(jj)*w(kk);
                        cont = cont+1;
                    end
                end
            end
            self.quadrature.Ng = Ngauss;
            self.quadrature.X = X;	% gauss integration points
            self.quadrature.W = W;	% gauss integration weights
            % INIZIALIZATION of some matrices (this should speedup
            % numerical integration)
            C = self.Material.get_stress_strain_matrix_3D;
            H = [1 0 0 0 0 0 0 0 0;
                0 0 0 0 1 0 0 0 0;
                0 0 0 0 0 0 0 0 1;
                0 1 0 1 0 0 0 0 0;
                0 0 1 0 0 0 1 0 0;
                0 0 0 0 0 1 0 1 0];
            self.initialization.A = zeros(6,9); % nonlinear strain
            self.initialization.G = zeros(9,60);% shape function derivatives
            self.initialization.Z = zeros(20);  % zero-matrix
            self.initialization.K = zeros(60);  % stiffness-element matrix
            self.initialization.F = zeros(60,1);% internal forces (element)
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
            Mel = zeros(60);
            for ii = 1:length(W)
                g = X(1,ii);
                h = X(2,ii);
                r = X(3,ii);
                we = W(ii); % weights
                N = self.shape_functions(g,h,r);
                [~,detJ] = shape_function_derivatives(self,g,h,r);
                NN(1,1:3:60) = N;
                NN(2,2:3:60) = N;
                NN(3,3:3:60) = N;
                % integration of K and M through GAUSS QUADRATURE
                Mel = Mel + (NN'*NN)*(we*detJ);
            end
            Mel = sparse(rho*Mel);
        end
        
        function [K,F] = tangent_stiffness_and_force(self,x)
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
                r = X(3,ii);
                we = W(ii); % weights
                [G,detJ,dH] = shape_function_derivatives(self,g,h,r);
                th  = G*displ;
                A = self.initialization.A;
                A(1,1)=th(1); A(4,1)=th(2); A(5,1)=th(3); A(2,2)=th(2); A(4,2)=th(1);
                A(6,2)=th(3); A(3,3)=th(3); A(5,3)=th(1); A(6,3)=th(2); A(1,4)=th(4);
                A(4,4)=th(5); A(5,4)=th(6); A(2,5)=th(5); A(4,5)=th(4); A(6,5)=th(6);
                A(3,6)=th(6); A(5,6)=th(4); A(6,6)=th(5); A(1,7)=th(7); A(4,7)=th(8);
                A(5,7)=th(9); A(2,8)=th(8); A(4,8)=th(7); A(6,8)=th(9); A(3,9)=th(9);
                A(5,9)=th(7); A(6,9)=th(8);
                % Green Strain tensor
                E = (H + 1/2*A)*th;
                % second Piola-Kirchhoff stress tensor
                s = C*E; % s = [S11 S22 S33 S12 S13 S23]
                S = [s(1) s(4) s(5); s(4) s(2) s(6); s(5) s(6) s(3)];
                Bnl = (H + A)*G;
                % functions to integrate over volume
                int_K1 = Bnl'*C*Bnl;
                HSH = dH'*S*dH;
                int_Ks = [HSH ZZ ZZ; ZZ HSH ZZ; ZZ ZZ HSH]; % (faster than blkdiag)
                int_K = (int_K1 + int_Ks)*detJ;
                int_F = (Bnl'*s)*detJ;
                % integration of K and F through Gauss quadrature
                K = K + we*int_K;
                F = F + we*int_F;
            end
        end
        
        function xe = extract_element_data(self,x)
            % x is a vector of full DOFs
            index = get_index(self.nodeIDs,self.nDOFPerNode);
            xe = x(index,:);
        end
        
        function F = get.uniformBodyForce(self)
            % _____________________________________________________________
            %
            % F = uniform_body_force(self,direction)
            % This function computes a load along direction=3(Z) by
            % dividing the load on the 20 nodes according to the element
            % volume (V/20) [it might not be the best way, but still...]
            %______________________________________________________________
            F = sparse(60,1);
            F(3:3:end) = self.vol/20; % uniformly distributed pressure on the structure
        end
        
        function [T2, globalSubs] = T2(self)
            % this function computes the 3-tensor corresponding to the 
            % quadratic component of the nonlinear internal force in 
            % global coordinates at the element level.
                        
            % global DOFs associated to the element nodes
            index = get_index(self.nodeIDs,self.nDOFPerNode);
            
            % location of each dimension of tensor in global DOFs
            globalSubs = {index, index, index};
                        
            X = self.quadrature.X;
            W = self.quadrature.W;

            C = self.initialization.C;  % constitutive law matrix
            H = self.initialization.H;  % Linear strain matrix: eps_l = H*th

            % Quadratic strain matrix: A = L.th, eps_quad = A*th
            L = tenzeros([6,9,9]);
            L(1,1,1)=1; L(4,2,1)=1; L(5,3,1)=1; 
            L(4,1,2)=1; L(2,2,2)=1; L(6,3,2)=1; 
            L(5,1,3)=1; L(6,2,3)=1; L(3,3,3)=1;
            L(1,4,4)=1; L(4,5,4)=1; L(5,6,4)=1; 
            L(4,4,5)=1; L(2,5,5)=1; L(6,6,5)=1; 
            L(5,4,6)=1; L(6,5,6)=1; L(3,6,6)=1;
            L(1,7,7)=1; L(4,8,7)=1; L(5,9,7)=1; 
            L(4,7,8)=1; L(2,8,8)=1; L(6,9,8)=1; 
            L(5,7,9)=1; L(6,8,9)=1; L(3,9,9)=1;

            m = self.nNodes*self.nDOFPerNode;
            Q3h = tenzeros([m,m,m]);
            for ii = 1:length(W)
                g = X(1,ii);
                h = X(2,ii);
                r = X(3,ii);
                we = W(ii); % weights
                [G,detJ] = shape_function_derivatives(self,g,h,r);

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
        
        function [T3, globalSubs] = T3(self)
            % this function computes the 4-tensor corresponding to the 
            % quadratic component of the nonlinear internal force in 
            % global coordinates at the element level.
                        
            % global DOFs associated to the element nodes
            index = get_index(self.nodeIDs,self.nDOFPerNode);
            
            % location of each dimension of tensor in global DOFs
            globalSubs = cell(4,1);
            globalSubs(:) = {index};
                        
            X = self.quadrature.X;
            W = self.quadrature.W;

            C = self.initialization.C;  % constitutive law matrix

            % Quadratic strain matrix: A = L.th, eps_quad = A*th
            L = tenzeros([6,9,9]);
            L(1,1,1)=1; L(4,2,1)=1; L(5,3,1)=1; 
            L(4,1,2)=1; L(2,2,2)=1; L(6,3,2)=1; 
            L(5,1,3)=1; L(6,2,3)=1; L(3,3,3)=1;
            L(1,4,4)=1; L(4,5,4)=1; L(5,6,4)=1; 
            L(4,4,5)=1; L(2,5,5)=1; L(6,6,5)=1; 
            L(5,4,6)=1; L(6,5,6)=1; L(3,6,6)=1;
            L(1,7,7)=1; L(4,8,7)=1; L(5,9,7)=1; 
            L(4,7,8)=1; L(2,8,8)=1; L(6,9,8)=1; 
            L(5,7,9)=1; L(6,8,9)=1; L(3,9,9)=1;

            m = self.nNodes*self.nDOFPerNode;
            T3 = tenzeros([m,m,m,m]);

            for ii = 1:length(W)
                g = X(1,ii);
                h = X(2,ii);
                r = X(3,ii);
                we = W(ii); % weights
                [G,detJ] = shape_function_derivatives(self,g,h,r);

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
            A = self.initialization.A;
            % Quadratic strain matrix: A = L.th, eps_quad = A*th
            L = zeros([3,4,4]);
            L(1,1,1)=1; L(4,2,1)=1; L(5,3,1)=1;
            L(4,1,2)=1; L(2,2,2)=1; L(6,3,2)=1;
            L(5,1,3)=1; L(6,2,3)=1; L(3,3,3)=1;
            L(1,4,4)=1; L(4,5,4)=1; L(5,6,4)=1;
            L(4,4,5)=1; L(2,5,5)=1; L(6,6,5)=1;
            L(5,4,6)=1; L(6,5,6)=1; L(3,6,6)=1;
            L(1,7,7)=1; L(4,8,7)=1; L(5,9,7)=1;
            L(4,7,8)=1; L(2,8,8)=1; L(6,9,8)=1;
            L(5,7,9)=1; L(6,8,9)=1; L(3,9,9)=1;
            for ii = 1:length(W)
                g = X(1,ii);
                h = X(2,ii);
                r = X(3,ii);
                we = W(ii); % weights
                [G,detJ] = shape_function_derivatives(self,g,h,r);
                th  = G*displ;
                A(1,1)=th(1); A(4,1)=th(2); A(5,1)=th(3); A(2,2)=th(2); A(4,2)=th(1);
                A(6,2)=th(3); A(3,3)=th(3); A(5,3)=th(1); A(6,3)=th(2); A(1,4)=th(4);
                A(4,4)=th(5); A(5,4)=th(6); A(2,5)=th(5); A(4,5)=th(4); A(6,5)=th(6);
                A(3,6)=th(6); A(5,6)=th(4); A(6,6)=th(5); A(1,7)=th(7); A(4,7)=th(8);
                A(5,7)=th(9); A(2,8)=th(8); A(4,8)=th(7); A(6,8)=th(9); A(3,9)=th(9);
                A(5,9)=th(7); A(6,9)=th(8);
                b1 = einsum('ijk,kl->ijl',L,G);
                b2 = einsum('ijk,il->jkl',b1,C*H*th);
                b = G'*b2;
                int_dK = G'*(H'*C*A + A'*C*H)*G + b;
                Kd = Kd + int_dK * detJ * we;
            end
        end
        
        % ANCILLARY FUNCTIONS _____________________________________________
        
        function V = get.vol(self)
            % volume is given by the integral of detJ (jacobian from 
            % isoparametric to physical space) over the volume of the
            % isoparametric element
            detJ = 0;
            W = self.quadrature.W;
            X = self.quadrature.X;
            for ii = 1 : length( w )
                g = X(1,ii);
                h = X(2,ii);
                r = X(3,ii);
                [~, detJ_i] = shape_function_derivatives(self,g,h,r);
                detJ = detJ + detJ_i * W(ii);
            end
            V = detJ;
        end
        
        function [G,detJ,dH] = shape_function_derivatives(self,g,h,r)
            %______________________________________________________________
            %
            % [G,detJ,dH] = G_HEX20(self,g,h,r)
            % G = shape function derivatives in physical coordinates, s.t.
            % th=G*p with th={ux uy uz vx vy vz wx wy wz}' (ux=du/dx...)
            % and p={u1,v1,w1,...,u20,v20,w20}' (nodal displacements).
            % detJ = det(J), J=jacobian
            %______________________________________________________________
            xyz = self.nodes;
            % shape functions derivatives (ABAQUS ORDER)
            dHn = [ ...
                (g/8 - 1/8)*(h - 1)*(r - 1) + ((h - 1)*(r - 1)*(g + h + r + 2))/8,   (g/8 + 1/8)*(h - 1)*(r - 1) - ((h - 1)*(r - 1)*(h - g + r + 2))/8, - (g/8 + 1/8)*(h + 1)*(r - 1) - ((h + 1)*(r - 1)*(g + h - r - 2))/8, - (g/8 - 1/8)*(h + 1)*(r - 1) - ((h + 1)*(r - 1)*(g - h + r + 2))/8, - (g/8 - 1/8)*(h - 1)*(r + 1) - ((h - 1)*(r + 1)*(g + h - r + 2))/8, - (g/8 + 1/8)*(h - 1)*(r + 1) - ((h - 1)*(r + 1)*(g - h + r - 2))/8, (g/8 + 1/8)*(h + 1)*(r + 1) + ((h + 1)*(r + 1)*(g + h + r - 2))/8, (g/8 - 1/8)*(h + 1)*(r + 1) + ((h + 1)*(r + 1)*(g - h - r + 2))/8, - (g/4 - 1/4)*(h - 1)*(r - 1) - ((g + 1)*(h - 1)*(r - 1))/4,                               (h/4 - 1/4)*(h + 1)*(r - 1), (g/4 - 1/4)*(h + 1)*(r - 1) + ((g + 1)*(h + 1)*(r - 1))/4,                                -(h/4 - 1/4)*(h + 1)*(r - 1), (g/4 - 1/4)*(h - 1)*(r + 1) + ((g + 1)*(h - 1)*(r + 1))/4,                                -(h/4 - 1/4)*(h + 1)*(r + 1), - (g/4 - 1/4)*(h + 1)*(r + 1) - ((g + 1)*(h + 1)*(r + 1))/4,                               (h/4 - 1/4)*(h + 1)*(r + 1),                                -(r/4 - 1/4)*(h - 1)*(r + 1),                               (r/4 - 1/4)*(h - 1)*(r + 1),                                -(r/4 - 1/4)*(h + 1)*(r + 1),                               (r/4 - 1/4)*(h + 1)*(r + 1);
                (g/8 - 1/8)*(h - 1)*(r - 1) + (g/8 - 1/8)*(r - 1)*(g + h + r + 2), - (g/8 + 1/8)*(h - 1)*(r - 1) - (g/8 + 1/8)*(r - 1)*(h - g + r + 2), - (g/8 + 1/8)*(h + 1)*(r - 1) - (g/8 + 1/8)*(r - 1)*(g + h - r - 2),   (g/8 - 1/8)*(h + 1)*(r - 1) - (g/8 - 1/8)*(r - 1)*(g - h + r + 2), - (g/8 - 1/8)*(h - 1)*(r + 1) - (g/8 - 1/8)*(r + 1)*(g + h - r + 2),   (g/8 + 1/8)*(h - 1)*(r + 1) - (g/8 + 1/8)*(r + 1)*(g - h + r - 2), (g/8 + 1/8)*(h + 1)*(r + 1) + (g/8 + 1/8)*(r + 1)*(g + h + r - 2), (g/8 - 1/8)*(r + 1)*(g - h - r + 2) - (g/8 - 1/8)*(h + 1)*(r + 1),                                -(g/4 - 1/4)*(g + 1)*(r - 1), (h/4 - 1/4)*(g + 1)*(r - 1) + ((g + 1)*(h + 1)*(r - 1))/4,                               (g/4 - 1/4)*(g + 1)*(r - 1), - (h/4 - 1/4)*(g - 1)*(r - 1) - ((g - 1)*(h + 1)*(r - 1))/4,                               (g/4 - 1/4)*(g + 1)*(r + 1), - (h/4 - 1/4)*(g + 1)*(r + 1) - ((g + 1)*(h + 1)*(r + 1))/4,                                -(g/4 - 1/4)*(g + 1)*(r + 1), (h/4 - 1/4)*(g - 1)*(r + 1) + ((g - 1)*(h + 1)*(r + 1))/4,                                -(r/4 - 1/4)*(g - 1)*(r + 1),                               (r/4 - 1/4)*(g + 1)*(r + 1),                                -(r/4 - 1/4)*(g + 1)*(r + 1),                               (r/4 - 1/4)*(g - 1)*(r + 1);
                (g/8 - 1/8)*(h - 1)*(r - 1) + (g/8 - 1/8)*(h - 1)*(g + h + r + 2), - (g/8 + 1/8)*(h - 1)*(r - 1) - (g/8 + 1/8)*(h - 1)*(h - g + r + 2),   (g/8 + 1/8)*(h + 1)*(r - 1) - (g/8 + 1/8)*(h + 1)*(g + h - r - 2), - (g/8 - 1/8)*(h + 1)*(r - 1) - (g/8 - 1/8)*(h + 1)*(g - h + r + 2),   (g/8 - 1/8)*(h - 1)*(r + 1) - (g/8 - 1/8)*(h - 1)*(g + h - r + 2), - (g/8 + 1/8)*(h - 1)*(r + 1) - (g/8 + 1/8)*(h - 1)*(g - h + r - 2), (g/8 + 1/8)*(h + 1)*(r + 1) + (g/8 + 1/8)*(h + 1)*(g + h + r - 2), (g/8 - 1/8)*(h + 1)*(g - h - r + 2) - (g/8 - 1/8)*(h + 1)*(r + 1),                                -(g/4 - 1/4)*(g + 1)*(h - 1),                               (h/4 - 1/4)*(g + 1)*(h + 1),                               (g/4 - 1/4)*(g + 1)*(h + 1),                                -(h/4 - 1/4)*(g - 1)*(h + 1),                               (g/4 - 1/4)*(g + 1)*(h - 1),                                -(h/4 - 1/4)*(g + 1)*(h + 1),                                -(g/4 - 1/4)*(g + 1)*(h + 1),                               (h/4 - 1/4)*(g - 1)*(h + 1), - (r/4 - 1/4)*(g - 1)*(h - 1) - ((g - 1)*(h - 1)*(r + 1))/4, (r/4 - 1/4)*(g + 1)*(h - 1) + ((g + 1)*(h - 1)*(r + 1))/4, - (r/4 - 1/4)*(g + 1)*(h + 1) - ((g + 1)*(h + 1)*(r + 1))/4, (r/4 - 1/4)*(g - 1)*(h + 1) + ((g - 1)*(h + 1)*(r + 1))/4];
            % jacobian
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
            dH = J1*dHn;        % derivatives in physical coordinates,
            % 3x20 matrix, [dNi_dx; dNi_dy; dNi_dz]
            % with i = 1...20
            G = self.initialization.G;
            G(1:3,1:3:60) = dH;
            G(4:6,2:3:60) = dH;
            G(7:9,3:3:60) = dH;
        end
        
    end % methods
    
    methods (Static)
        
        function N = shape_functions(g,h,r)
            % N = shape_functions(g,h,r)
            % SHAPE FUNCTIONS FOR A 20-NODED HEXAHEDRON
            % - see Abaqus documentation:
            % Abaqus theory guide > Elements > Continuum elements > ...
            % ... > 3.2.4 Solid isoparametric quadrilaterals and hexahedra
            N = [ -1/8*(1-g)*(1-h)*(1-r)*(2+g+h+r); % n1
                -1/8*(1+g)*(1-h)*(1-r)*(2-g+h+r);   % n2
                -1/8*(1+g)*(1+h)*(1-r)*(2-g-h+r);   % n3
                -1/8*(1-g)*(1+h)*(1-r)*(2+g-h+r);   % n4
                -1/8*(1-g)*(1-h)*(1+r)*(2+g+h-r);   % n5
                -1/8*(1+g)*(1-h)*(1+r)*(2-g+h-r);   % n6
                -1/8*(1+g)*(1+h)*(1+r)*(2-g-h-r);   % n7
                -1/8*(1-g)*(1+h)*(1+r)*(2+g-h-r);   % n8
                +1/4*(1-g)*(1+g)*(1-h)*(1-r);       % n9
                +1/4*(1-h)*(1+h)*(1+g)*(1-r);       % n10
                +1/4*(1-g)*(1+g)*(1+h)*(1-r);       % n11
                +1/4*(1-h)*(1+h)*(1-g)*(1-r);       % n12
                +1/4*(1-g)*(1+g)*(1-h)*(1+r);       % n13
                +1/4*(1-h)*(1+h)*(1+g)*(1+r);       % n14
                +1/4*(1-g)*(1+g)*(1+h)*(1+r);       % n15
                +1/4*(1-h)*(1+h)*(1-g)*(1+r);       % n16
                +1/4*(1-r)*(1+r)*(1-g)*(1-h);       % n17
                +1/4*(1-r)*(1+r)*(1+g)*(1-h);       % n18
                +1/4*(1-r)*(1+r)*(1+g)*(1+h);       % n19
                +1/4*(1-r)*(1+r)*(1-g)*(1+h)];      % n20
        end
        
        function X = natural_coordinates
            X = [ ...
                -1 -1 -1
                +1 -1 -1
                +1 +1 -1
                -1 +1 -1
                -1 -1 +1
                +1 -1 +1
                +1 +1 +1
                -1 +1 +1
                 0 -1 -1
                +1  0 -1
                 0 +1 -1
                -1  0 -1
                 0 -1 +1
                +1  0 +1
                 0 +1 +1
                -1  0 +1
                -1 -1  0
                +1 -1  0
                +1 +1  0
                -1 +1  0];
        end
        
    end
    
end % classdef


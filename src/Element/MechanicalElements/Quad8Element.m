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
            [x, w] = lgwt(Ngauss,-1,1);
            self.quadrature.Ng = Ngauss;
            self.quadrature.X = x;      % gauss integration points
            self.quadrature.W = w;      % gauss integration weights
            self.thickness = thickness;
            % INIZIALIZATION of some matrices (this should speedup
            % numerical integration)
            C = self.Material.get_stress_strain_matrix_2D;
            H = [1 0 0 0;
                0 0 0 1;
                0 1 1 0];
            % Quadratic strain matrix: A = L.th, eps_quad = A*th
            L(1,1,1)=1; L(3,2,1)=1; L(3,1,2)=1; L(2,2,2)=1; L(1,2,3)=1;
            L(3,4,3)=1; L(3,2,4)=1; L(2,4,4)=1;
            L = tensor(L);
            
            self.initialization.A = zeros(3,4); % nonlinear strain
            self.initialization.G = zeros(4,16);% shape function derivatives
            self.initialization.Z = zeros(8);   % zero-matrix
            self.initialization.K = zeros(16);  % stiffness-element matrix
            self.initialization.F = zeros(16,1);% internal forces (element)
            self.initialization.C = C;          % constitutive law matrix
            self.initialization.H = H;          % linear strain
            self.initialization.L=L;            %Quadratic strain matrix: A = L.th,
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
            Mel = zeros(16);
            for ii = 1:self.quadrature.Ng
                for jj = 1:self.quadrature.Ng
                    g = X(ii);
                    h = X(jj);
                    % shape functions and detJ (ABAQUS ORDER)
                    N = self.shape_functions(g,h);
                    [~,detJ] = shape_function_derivatives(self,g,h);
                    NN = kron(N', eye(2));
                    % integration of K and M through GAUSS QUADRATURE
                    Mel = Mel + ( NN' * NN )*( W(ii) * W(jj) * detJ );
                end
            end
            Mel = sparse(rho*Mel);
        end
        function Cel = damping_matrix(self,alfa,beta,u0)
            [Kel Fel]=self.tangent_stiffness_and_force(u0);
            
            Cel=alfa*self.mass_matrix()+beta*Kel;
        end
        function Kel = stiffness_matrix(self,u0)
            [Kel Fel]=self.tangent_stiffness_and_force(u0);
            
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
            for ii = 1:self.quadrature.Ng
                for jj = 1:self.quadrature.Ng
                    g = X(ii);
                    h = X(jj);
                    we = W(ii)*W(jj); % weights
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
        end
        
        function xe = extract_element_data(self, x)
            % x is a vector of full DOFs
            index = get_index(self.nodeIDs, self.nDOFPerNode);
            xe = x(index,:);
        end
        
        function f = uniform_body_force(self,~)
            % _____________________________________________________________
            %
            % F = uniform_body_force(self,direction)
            % This function computes a load along direction=2(Y) by
            % dividing the load on the 8 nodes according to the element
            % area (A/8) [it might not be the best way, but still...]
            %______________________________________________________________
            f = sparse(16,1);
            f(2:2:end) = self.area/8; % uniformly distributed pressure on the structure
        end
        
        % ANCILLARY FUNCTIONS _____________________________________________
        
        function A = get_area(self)
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
        
        function dKdq= stiffness_derivative(self,u0,VMs)
            
            X=self.quadrature.X;  %g=X(ii) h=X(jj);
            %             g  = X(1);
            %             h = X(2);
            W=self.quadrature.W;  %weights
            H=self.initialization.H;
            C=self.Material.get_stress_strain_matrix_2D;
            A=self.initialization.A;
            
            K = self.stiffness_matrix(u0);
            M=self.mass_matrix();
            
            %           %nnodes = size(xyz,1); %get Phi_i_elem
            %             [Phi,D] = eigs(K,M,Nm,'SM'); [f0,ind] =
            %             sort(sqrt(diag(D))/2/pi); VMs = Phi(:,ind);
            
            %normalize
            %             for ii = 1:Nm
            %                 VMs(:,ii) =
            %                 VMs(:,ii)/max(sqrt(sum(VMs(:,ii).^2,2)));
            %             end
            Phi_i_el=VMs(self.iDOFs);
            
            %             [G,detJ,dH] = self.shape_function_derivatives(g,h);
            dKdq=0;
            for ii = 1:self.quadrature.Ng
                for jj = 1:self.quadrature.Ng
                    g = X(ii);
                    h = X(jj);
                    we = W(ii)*W(jj); % weights
                    [G,detJ,dH] = shape_function_derivatives(self,g,h);
                    
                    
                    G=sparse(G);
                    % dK(q)/dq_i
                    th = G*Phi_i_el; % displacement derivatives vector
                    
                    A =	[th(1)	0     th(3) 0;
                        0    	th(2) 0     th(4);
                        th(2)	th(1) th(4) th(3)];
                    int_dKdq = G'*(H'*C*A + 2*A'*C*H)*G*detJ;
                    dKdq=dKdq+we*int_dKdq;
                    
                end
            end
        end
        
        function [Q3 globalSubs]= tensor_Q3(self)
            X=self.quadrature.X;  %g=X(ii) h=X(jj);
            W=self.quadrature.W;  %weights
            H=self.initialization.H;
            C=self.Material.get_stress_strain_matrix_2D;
            %             A=self.initialization.A;
            L=self.initialization.L;
            
                         % Phi_i_el=V(self.iDOFs);
            
            m=numel(self.nodes);
            Q3h = zeros(m,m,m);
            
            for ii = 1:self.quadrature.Ng
                for jj = 1:self.quadrature.Ng
                    g = X(ii);
                    h = X(jj);
                    we = W(ii)*W(jj); % weights
                    [G,detJ,dH] = shape_function_derivatives(self,g,h);
                    
                    %                     GHC = tensor((C*H*G)');
                    %
                    %                     C = tensor(C);
                    %                     G = tensor(G);
                    %
                    %                     LGG = ttt(ttt(L,G,3,1),G,2,1);
                    %
                    %                     Q3h_int = ttt(GHC,LGG,2,1);
                    %
                    Q3h_int = self.tensors_Q3_hat(G,H,L,C);
                    
                    
                    Q3h = Q3h + Q3h_int*detJ*we;
                    
                    
                end
            end
            % build third order tensors using Q3h
            Q3ht = permute(Q3h,[3 2 1]);
            Q3 = Q3h./2 + Q3ht;
            
%                            Vt = tensor(Phi_i_el');    % reduction basis (transpose)
%                            V  = tensor(Phi_i_el);     % reduction basis
% %             
%             % third order tensor projections
                          %  Q3  = ttt(ttt(ttt(Vt,Q3 ,2,1),V,3,1),V ,2,1);
%             %               % Taken from rectangular plate file
            % global DOFs associated to the element nodes
            index = get_index(self.nodeIDs,self.nDOFPerNode);
            
            % location of each dimension of tensor in global DOFs
            globalSubs = {index, index, index};
            
            
            
        end
        
        function [Q3h] = tensors_Q3_hat(self,G,H,L,C)
            
            GHC = tensor((C*H*G)');
            C = tensor(C);
            G = tensor(G);
            L=tensor(L);
            LGG = tensor(ttt(ttt(L,G,3,1),G,2,1));
            
            Q3h = ttt(GHC,LGG,2,1);
        end
        function [Q4 globalSubs]= tensor_Q4(self)
            X=self.quadrature.X;  %g=X(ii) h=X(jj);
            W=self.quadrature.W;  %weights
            H=self.initialization.H;
            C=self.Material.get_stress_strain_matrix_2D;
            %             A=self.initialization.A;
            L=self.initialization.L;
            
            %Phi_i_el=V(self.iDOFs);
            
            m=numel(self.nodes);
            Q4h = zeros(m,m,m,m);
            for ii = 1:self.quadrature.Ng
                for jj = 1:self.quadrature.Ng
                    g = X(ii);
                    h = X(jj);
                    we = W(ii)*W(jj); % weights
                    [G,detJ,dH] = shape_function_derivatives(self,g,h);
                    %                     G=sparse(G);
                    %
                    %                     GHC = tensor((C*H*G)');
                    %                     C = tensor(C);
                    %                     G = tensor(G);
                    %
                    %                     LGG = ttt(ttt(L,G,3,1),G,2,1);
                    %
                    %
                    %                     Q4h_int = ttt(ttt(permute(LGG,[2 1 3]),C,2,1),LGG,3,1);
                    %
                    Q4h_int = self.tensors_Q4_hat(G,H,L,C);
                    
                    
                    Q4h = Q4h + Q4h_int*detJ*we;
                    
                end
            end
            Q4=Q4h./2;
%             Vt = tensor(V');    % reduction basis (transpose)
%             V  = tensor(V);     % reduction basis
%             
%             
%             % fourth order tensor projections
%             Q4   = ttt(ttt(ttt(ttt(Vt,Q4  ,2,1),V,4,1),V ,3,1),V ,2,1);
%             
%             
            % global DOFs associated to the element nodes
            index = get_index(self.nodeIDs,self.nDOFPerNode);
            
            % location of each dimension of tensor in global DOFs
            globalSubs = cell(4,1);
            globalSubs(:) = {index};
        end
        function [Q4h] = tensors_Q4_hat(self,G,H,L,C)
            
           
            C = tensor(C);
            G = tensor(G);
            L=tensor(L);

            LGG = tensor(ttt(ttt(L,G,3,1),G,2,1));
            
            
            Q4h = ttt(ttt(permute(LGG,[2 1 3]),C,2,1),LGG,3,1);
            
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


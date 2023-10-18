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
                      	% 2x6 matrix, [dNi_dx; dNi_dy]
                       	% with i = 1...6
            G = self.initialization.G;
            G(1:2,1:2:end) = dH;
            G(3:4,2:2:end) = dH;
        end
        
        % HYDRODYNAMIC FUNCTIONS __________________________________________

        % Drag force tensors
        function T = Te1(self, specificFace, rho)
            % _____________________________________________________________
            % Returns one of the 2nd order tensors (i.e., a matrix) stemming from 
            % the drag force. These forces are acting on the
            % ``specificFace'' of an element. The forces are divided by 2 
            % and applied to the nodes at of the considered face. 
            % _____________________________________________________________

            if specificFace(1) == 1
                T = T1_conf1(rho,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
            elseif specificFace(1) == 2
                T = T1_conf2(rho,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
            else
                T = T1_conf3(rho,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
            end
            
            if specificFace(2) ~= 0                               
                if specificFace(2) == 1
                    T = T1_conf1(rho,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
                elseif specificFace(2) == 2
                    T = T1_conf2(rho,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
                else
                    T = T1_conf3(rho,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
                end
            end
        end

        function T = Teu2(self, specificFace, rho)
            % _____________________________________________________________
            % Returns one of the 2nd order tensors (i.e., a matrix) stemming from 
            % the drag force. These forces are acting on the
            % ``specificFace'' of an element. The forces are divided by 2 
            % and applied to the nodes at of the considered face. 
            % _____________________________________________________________

            if specificFace(1) == 1
                T = Tu2_conf1(rho,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
            elseif specificFace(1) == 2
                T = Tu2_conf2(rho,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
            else
                T = Tu2_conf3(rho,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
            end
            
            if specificFace(2) ~= 0                               
                if specificFace(2) == 1
                    T = Tu2_conf1(rho,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
                elseif specificFace(2) == 2
                    T = Tu2_conf2(rho,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
                else
                    T = Tu2_conf3(rho,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
                end
            end
        end

        function T = Teu3(self, specificFace, vwater, rho, c)
            % _____________________________________________________________
            % Returns one of the 3rd order tensors stemming from 
            % the drag and the thrust force. These forces are acting on the
            % ``specificFace'' of an element. The forces are divided by 2 
            % and applied to the nodes at of the considered face. 
            % _____________________________________________________________
            [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(1));
            inv = normal_vector_inversion(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));

            if specificFace(1) == 1
                T = Tu3_conf1_V4(rho,c,vwater(1),vwater(2),inv,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
            elseif specificFace(1) == 2
                T = Tu3_conf2_V4(rho,c,vwater(1),vwater(2),inv,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
            else
                T = Tu3_conf3_V4(rho,c,vwater(1),vwater(2),inv,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
            end
            
            if specificFace(2) ~= 0
                [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(2));
                inv = normal_vector_inversion(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
                
                if specificFace(2) == 1
                    T = Tu3_conf1_V4(rho,c,vwater(1),vwater(2),inv,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
                elseif specificFace(2) == 2
                    T = Tu3_conf2_V4(rho,c,vwater(1),vwater(2),inv,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
                else
                    T = Tu3_conf3_V4(rho,c,vwater(1),vwater(2),inv,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
                end
            end

            T = permute(T,[3 1 2]);

        end

        function T = Teudot2(self, specificFace, rho)
            % _____________________________________________________________
            % Returns one of the 2nd order tensors (i.e., a matrix) stemming from 
            % the drag force. These forces are acting on the
            % ``specificFace'' of an element. The forces are divided by 2 
            % and applied to the nodes at of the considered face. 
            % _____________________________________________________________
           
            if specificFace(1) == 1
                T = Tudot2_conf1(rho,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
            elseif specificFace(1) == 2
                T = Tudot2_conf2(rho,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
            else
                T = Tudot2_conf3(rho,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
            end
            
            if specificFace(2) ~= 0
               
                if specificFace(2) == 1
                    T = Tudot2_conf1(rho,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
                elseif specificFace(2) == 2
                    T = Tudot2_conf2(rho,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
                else
                    T = Tudot2_conf3(rho,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
                end
            end

        end

        function T = Teudot3(self, specificFace, vwater, rho, c)
            % _____________________________________________________________
            % Returns one of the 2nd order tensors (i.e., a matrix) stemming from 
            % the drag and the thrust force. These forces are acting on the
            % ``specificFace'' of an element. The forces are divided by 2 
            % and applied to the nodes at of the considered face. 
            % _____________________________________________________________
            [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(1));
            inv = normal_vector_inversion(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));

            if specificFace(1) == 1
                T = Tudot3_conf1_V4(rho,c,vwater(1),vwater(2),inv,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
            elseif specificFace(1) == 2
                T = Tudot3_conf2_V4(rho,c,vwater(1),vwater(2),inv,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
            else
                T = Tudot3_conf3_V4(rho,c,vwater(1),vwater(2),inv,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
            end
            
            if specificFace(2) ~= 0
                [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(2));
                inv = normal_vector_inversion(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
                
                if specificFace(2) == 1
                    T = Tudot3_conf1_V4(rho,c,vwater(1),vwater(2),inv,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
                elseif specificFace(2) == 2
                    T = Tudot3_conf2_V4(rho,c,vwater(1),vwater(2),inv,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
                else
                    T = Tudot3_conf3_V4(rho,c,vwater(1),vwater(2),inv,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
                end
            end

        end
        
        function T = Teuu3(self, specificFace, rho)
            % _____________________________________________________________
            % Returns one of the 3rd order tensors stemming from 
            % the drag force. These forces are acting on the
            % ``specificFace'' of an element. The forces are divided by 2 
            % and applied to the nodes at of the considered face. 
            % _____________________________________________________________
                   
            if specificFace(1) == 1
                T = Tuu3_conf1(rho,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
            elseif specificFace(1) == 2
                T = Tuu3_conf2(rho,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
            else
                T = Tuu3_conf3(rho,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
            end
            
            if specificFace(2) ~= 0
                                
                if specificFace(2) == 1
                    T = Tuu3_conf1(rho,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
                elseif specificFace(2) == 2
                    T = Tuu3_conf2(rho,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
                else
                    T = Tuu3_conf3(rho,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
                end
            end

            T = permute(T,[3 1 2]);
        end

        function T = Teuu4(self, specificFace, vwater, rho, c)
            % _____________________________________________________________
            % Returns one of the 4th order tensors stemming from 
            % the drag and the thrust force. These forces are acting on the
            % ``specificFace'' of an element. The forces are divided by 2 
            % and applied to the nodes at of the considered face. 
            % _____________________________________________________________
            [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(1));
            inv = normal_vector_inversion(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
            
            if specificFace(1) == 1
                T = Tuu4_conf1_V4(rho,c,vwater(1),vwater(2),inv,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));      
            elseif specificFace(1) == 2
                T = Tuu4_conf2_V4(rho,c,vwater(1),vwater(2),inv,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
            else
                T = Tuu4_conf3_V4(rho,c,vwater(1),vwater(2),inv,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
            end
            
            if specificFace(2) ~= 0
                [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(2));
                inv = normal_vector_inversion(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
                
                if specificFace(2) == 1
                    T = Tuu4_conf1_V4(rho,c,vwater(1),vwater(2),inv,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));      
                elseif specificFace(2) == 2
                    T = Tuu4_conf2_V4(rho,c,vwater(1),vwater(2),inv,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
                else
                    T = Tuu4_conf3_V4(rho,c,vwater(1),vwater(2),inv,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
                end
            end
    
            T = permute(T,[4 3 1 2]);
        end

        function T = Teuudot3(self, specificFace, rho)
            % _____________________________________________________________
            % Returns one of the 3rd order tensors stemming from 
            % the drag force. These forces are acting on the
            % ``specificFace'' of an element. The forces are divided by 2 
            % and applied to the nodes at of the considered face. 
            % _____________________________________________________________
            
            if specificFace(1) == 1
                T = Tuudot3_conf1(rho,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
            elseif specificFace(1) == 2
                T = Tuudot3_conf2(rho,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
            else
                T = Tuudot3_conf3(rho,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
            end
            
            if specificFace(2) ~= 0
                
                if specificFace(2) == 1
                    T = Tuudot3_conf1(rho,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
                elseif specificFace(2) == 2
                   T = Tuudot3_conf2(rho,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
                else
                   T = Tuudot3_conf3(rho,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
                end
            end

            T = permute(T,[3 1 2]);
        end

        function T = Teuudot4(self, specificFace, vwater, rho, c)
            % _____________________________________________________________
            % Returns one of the 4th order tensors stemming from 
            % the drag and the thrust force. These forces are acting on the
            % ``specificFace'' of an element. The forces are divided by 2 
            % and applied to the nodes at of the considered face. 
            % _____________________________________________________________
            [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(1));
            inv = normal_vector_inversion(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
            if specificFace(1) == 1
                T = Tuudot4_conf1_V4(rho,c,vwater(1),vwater(2),inv,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));      
            elseif specificFace(1) == 2
                T = Tuudot4_conf2_V4(rho,c,vwater(1),vwater(2),inv,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
            else
                T = Tuudot4_conf3_V4(rho,c,vwater(1),vwater(2),inv,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
            end
            
            if specificFace(2) ~= 0
                [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(2));
                inv = normal_vector_inversion(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
                
                if specificFace(2) == 1
                    T = Tuudot4_conf1_V4(rho,c,vwater(1),vwater(2),inv,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));      
                elseif specificFace(2) == 2
                    T = Tuudot4_conf2_V4(rho,c,vwater(1),vwater(2),inv,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
                else
                    T = Tuudot4_conf3_V4(rho,c,vwater(1),vwater(2),inv,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
                end
            end

            T = permute(T,[4 3 1 2]);
        end

        function T = Teudotudot3(self, specificFace, rho)
            % _____________________________________________________________
            % Returns one of the 3rd order tensors stemming from 
            % the drag force. These forces are acting on the
            % ``specificFace'' of an element. The forces are divided by 2 
            % and applied to the nodes at of the considered face. 
            % _____________________________________________________________
                        
            if specificFace(1) == 1
                T = Tudotudot3_conf1(rho,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
            elseif specificFace(1) == 2
                T = Tudotudot3_conf2(rho,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
            else
                T = Tudotudot3_conf3(rho,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
            end
            
            if specificFace(2) ~= 0
                                
                if specificFace(2) == 1
                    T = Tudotudot3_conf1(rho,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
                elseif specificFace(2) == 2
                    T = Tudotudot3_conf2(rho,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
                else
                    T = Tudotudot3_conf3(rho,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));
                end
            end

            T = permute(T,[3 1 2]);
        end

        function T = Teudotudot4(self, specificFace, vwater, rho, c)
            % _____________________________________________________________
            % Returns one of the 3rd order tensors stemming from 
            % the drag and the thrust force. These forces are acting on the
            % ``specificFace'' of an element. The forces are divided by 2 
            % and applied to the nodes at of the considered face. 
            % _____________________________________________________________
            [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(1));
            inv = normal_vector_inversion(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
            if specificFace(1) == 1
                T = Tudotudot4_conf1_V4(rho,c,vwater(1),vwater(2),inv,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));      
            elseif specificFace(1) == 2
                T = Tudotudot4_conf2_V4(rho,c,vwater(1),vwater(2),inv,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
            else
                T = Tudotudot4_conf3_V4(rho,c,vwater(1),vwater(2),inv,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
            end
            
            if specificFace(2) ~= 0
                [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(2));
                inv = normal_vector_inversion(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
                
                if specificFace(2) == 1
                    T = Tudotudot4_conf1_V4(rho,c,vwater(1),vwater(2),inv,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));      
                elseif specificFace(2) == 2
                    T = Tudotudot4_conf2_V4(rho,c,vwater(1),vwater(2),inv,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
                else
                    T = Tudotudot4_conf3_V4(rho,c,vwater(1),vwater(2),inv,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
                end
            end

            T = permute(T,[4 3 1 2]);
        end
        
        % Early stage and test functions
        function F = drag_force_full(self, specificFace, vwater, rho, c, u, ud)
            % _____________________________________________________________
            % Returns the drag force F_drag acting the face
            % ``specificFace'' of an element. The force is divided by 2 and
            % applied to the nodes at of the considered face.
            % _____________________________________________________________
            F = sparse(self.nelDOFs,1);
            [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(1));

            % Compute drag force
            Fdrag = compute_drag_force(startNode, endNode, nextNode);
            % Apply drag force on nodes
            F = apply_force(self, F, Fdrag, startNode, endNode);

            if specificFace(2) ~= 0
                [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(2));

                % Compute drag force
                Fdrag = compute_drag_force(startNode, endNode, nextNode);
                % Apply drag force on nodes
                F = apply_force(self, F, Fdrag, startNode, endNode);
            end

            function Fdrag = compute_drag_force(startNode, endNode, nextNode)
                % _________________________________________________________
                % Computes the drag force vector for a specific surface
                %__________________________________________________________
                % Terms needed for the drag force - REPLACE VREL (to do) !!
                n = normal_vector(self, self.nodes(startNode,:)+u(startNode,:), self.nodes(endNode,:)+u(endNode,:), self.nodes(nextNode,:)+u(nextNode,:));  
                length = norm(self.nodes(endNode,:)+u(endNode,:)-self.nodes(startNode,:)-u(startNode,:));
                A = length; % in 2D, area proportional force is length proportional
                ns = 2; % 2D surfaces of TRI3 are composed of 2 elements
                vrel = vwater - 1/ns*(ud(startNode,:)+ud(endNode,:)); 
                d = vrel/norm(vrel);
                Cd = 2*(n.'*d)^2;
                
                % Drag force
                Fdrag = 1/2*rho*A*Cd*norm(vrel)^2*d;
            end 
            
        end
 
        function F = force_length_prop_skin_normal(self,specificFace)
            % _____________________________________________________________
            % Applies a force proportional to the the length of a skin
            % element on the 2 corresponding nodes, in the direction normal
            % to the face considered. Robust to elements having up to 2 
            % skin faces.
            % _____________________________________________________________
            F = sparse(self.nelDOFs,1);
            [startNode, endNode, nextNode] = get_node_from_face(self, specificFace(1));
            n = normal_vector(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));

            length = norm(self.nodes(endNode,:)-self.nodes(startNode,:));
            Force = length*n;
            F = apply_force(self, F, Force, startNode, endNode);

            if specificFace(2) ~= 0
                [startNode, endNode, nextNode] = get_node_from_face(self, specificFace(2));
                n = normal_vector(self,self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));

                length = norm(self.nodes(endNode,:)-self.nodes(startNode,:));
                Force = length*n;
                F = apply_force(self, F, Force, startNode, endNode);
            end
        end 
        
        function F = apply_force(~, Fbase, FtoApply, startNode, endNode)
            % _____________________________________________________________
            % Returns the addition of the force to apply FtoApply and the
            % base force Fbase.
            % _____________________________________________________________
            F = Fbase;
            F(startNode*2-1) = Fbase(startNode*2-1) + FtoApply(1)/2;    % x-coordinate
            F(startNode*2) = Fbase(startNode*2) + FtoApply(2)/2;        % y-coordinate
            F(endNode*2-1) = Fbase(endNode*2-1) + FtoApply(1)/2;        % x-coordinate
            F(endNode*2) = Fbase(endNode*2) + FtoApply(2)/2;            % y-coordinate
        end 
               
        function F = drag_force(self, specificFace, vwater, rho)
            % _____________________________________________________________
            % Returns the drag force F_drag acting the face
            % ``specificFace'' of an element. The force is divided by 2 and
            % applied to the nodes at of the considered face. Use of 1st
            % order tensors.
            % _____________________________________________________________
            [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(1));
            length = norm(self.nodes(endNode,:)-self.nodes(startNode,:));
            [r11,r12,r21,r22] = normal_vec_rot_matrix(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
            udot = [1 0 1 0 1 0].';
           
            if specificFace(1) == 1
                F = Te_conf_1_drag(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).' ...
                    + Tudote_conf_1_drag(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0)*udot;
            elseif specificFace(1) == 2
                F = Te_conf_2_drag(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).'...
                    + Tudote_conf_2_drag(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0)*udot;
            else
                F = Te_conf_3_drag(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).'...
                    + Tudote_conf_3_drag(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0)*udot;
            end
             %disp(F)
            if specificFace(2) ~= 0
                [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(2));
                length = norm(self.nodes(endNode,:)-self.nodes(startNode,:));
                [r11,r12,r21,r22] = normal_vec_rot_matrix(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
                if specificFace(2) == 1
                    F = F + Te_conf_1_drag(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).'...
                        + Tudote_conf_1_drag(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0)*udot;
                elseif specificFace(2) == 2
                    F = F + Te_conf_2_drag(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).'...
                        + Tudote_conf_2_drag(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0)*udot;
                else
                    F = F + Te_conf_3_drag(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).'...
                        + Tudote_conf_3_drag(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0)*udot;
                end
            end
            
            
        end
        
        % Reactive force (thrust force) functions
        function F = tail_pressure_force(self, tailNodeIndexInElement, normalisation, mTilde, u, ud)
            % _____________________________________________________________ 
            % Compute the pressure force at the tail f=0.5*m*v_perp^2*n
            % _____________________________________________________________
            F = zeros(self.nelDOFs, 1);
            
            % get node deformation u and its velocity ud (size 3 x 2)
            u2Columns = Tri3Element.vector_to_matrix(u); 
            ud2Columns = Tri3Element.vector_to_matrix(ud);
            uElement = u2Columns(self.nodeIDs,:);
            udElement = ud2Columns(self.nodeIDs,:);
            nodesPos = self.nodes + uElement; % position = base node position + displacement
            
            % convert back to 6 x 1 (x y x y x y)
            nodesPos = reshape(nodesPos',[],1);  
            uElement = reshape(uElement.',[],1);
            udElement = reshape(udElement',[],1);

            % get matrix corresponding to the configuration
            if tailNodeIndexInElement(1)==1         % conf 1-2 and 1-3
                A = A_conf1;
                if tailNodeIndexInElement(2)==2
                    B = B_conf1_2;
                else
                    B = B_conf1_3;
                end
            elseif tailNodeIndexInElement(1)==2     % conf 2-1 and 2-3
                A = A_conf2;
                if tailNodeIndexInElement(2)==1
                    B = B_conf2_1;
                else
                    B = B_conf2_3;
                end
            else                                    % conf 3-1 and 3-2
                A = A_conf3;
                if tailNodeIndexInElement(2)==1
                    B = B_conf3_1;
                else
                    B = B_conf3_2;
                end
            end
            
            % compute and apply force to tail node
            R = [0 -1 0 0 0 0;
                 1 0 0 0 0 0;
                 0 0 0 -1 0 0;
                 0 0 1 0 0 0;
                 0 0 0 0 0 1;
                 0 0 0 0 -1 0];     % 90 degrees rotation counterclock-wise
            t = B*nodesPos;
            n = R*t;
            v_perp = dot(A*udElement,n);
            F = 0.5*mTilde*v_perp.^2*normalisation*t;
        end

        function F = spine_momentum_force(self, spineNodeIndexInElement, normalisation, mTilde, u, ud, udd)
            % _____________________________________________________________ 
            % Compute the pressure force at the tail f=0.5*m*v_perp^2*n
            % _____________________________________________________________
            F = zeros(self.nelDOFs, 1);

            % get node deformation u and its velocity ud (size 3 x 2)
            u2Columns = Tri3Element.vector_to_matrix(u); 
            ud2Columns = Tri3Element.vector_to_matrix(ud);
            udd2Columns = Tri3Element.vector_to_matrix(udd);
            uElement = u2Columns(self.nodeIDs,:);
            udElement = ud2Columns(self.nodeIDs,:);
            uddElement = udd2Columns(self.nodeIDs,:);
            nodesPos = self.nodes + uElement; % position = base node position + displacement

            % convert back to 6 x 1 (x y x y x y)
            nodesPos = reshape(nodesPos',[],1);  
            uElement = reshape(uElement.',[],1);
            udElement = reshape(udElement',[],1);
            uddElement = reshape(uddElement',[],1);

            % get matrix corresponding to the configuration
            if spineNodeIndexInElement(1)==1         % conf 1-2 and 1-3
                A = A_conf1;
                if spineNodeIndexInElement(2)==2
                    B = B_conf1_2;
                else
                    B = B_conf1_3;
                end
            elseif spineNodeIndexInElement(1)==2     % conf 2-1 and 2-3
                A = A_conf2;
                if spineNodeIndexInElement(2)==1
                    B = B_conf2_1;
                else
                    B = B_conf2_3;
                end
            else                                    % conf 3-1 and 3-2
                A = A_conf3;
                if spineNodeIndexInElement(2)==1
                    B = B_conf3_1;
                else
                    B = B_conf3_2;
                end
            end

            % compute and apply force to tail node
            R = [0 -1 0 0 0 0;
                 1 0 0 0 0 0;
                 0 0 0 -1 0 0;
                 0 0 1 0 0 0;
                 0 0 0 0 0 1;
                 0 0 0 0 -1 0];     % 90 degrees rotation counterclock-wise

            F = ((A*uddElement).'*R*B*nodesPos + (A*udElement).'*R*B*udElement)*R*B*nodesPos + (A*udElement).'*R*B*nodesPos*R*B*udElement;
            % F = ((A*uddElement).'*R*B*nodesPos)*R*B*nodesPos ;

            F = - mTilde*normalisation*F;

            T = tensor(-mTilde*normalisation*einsum('mJ,mi,iK,Is,sL->IJKL',A,R,B,R,B));

            F = double(ttv(ttv(ttv(T,uddElement,2),nodesPos,2),nodesPos,2) + ttv(ttv(ttv(T,udElement,2),udElement,2),nodesPos,2)+ ttv(ttv(ttv(T,udElement,2),nodesPos,2),udElement,2));
        end

        function T = spine_momentum_tensor(self, tailNodeIndexInElement, normalisation, mTilde)
            % _____________________________________________________________ 
            % Compute the spine momentum change tensor needed to compute 
            % d/dt \int_0^l mTilde(a) v_\perp n da
            % _____________________________________________________________

            % get matrix corresponding to the configuration
            if tailNodeIndexInElement(1)==1         % conf 1-2 and 1-3
                A = A_conf1;
                if tailNodeIndexInElement(2)==2
                    B = B_conf1_2;
                else
                    B = B_conf1_3;
                end
            elseif tailNodeIndexInElement(1)==2     % conf 2-1 and 2-3
                A = A_conf2;
                if tailNodeIndexInElement(2)==1
                    B = B_conf2_1;
                else
                    B = B_conf2_3;
                end
            else                                    % conf 3-1 and 3-2
                A = A_conf3;
                if tailNodeIndexInElement(2)==1
                    B = B_conf3_1;
                else
                    B = B_conf3_2;
                end
            end

            % compute tensor
            R = [0 -1 0 0 0 0;
                 1 0 0 0 0 0;
                 0 0 0 -1 0 0;
                 0 0 1 0 0 0;
                 0 0 0 0 0 1;
                 0 0 0 0 -1 0];     % 90 degrees rotation counterclock-wise

            T = -mTilde*normalisation*einsum('mJ,mi,iK,Is,sL->IJKL',A,R,B,R,B);


            % % test
            % u = reshape(nodes.',[],1);
            % ud = zeros(size(u));
            % ud(2:2:end) = 4;
            % ud(1:2:end) = 0;
            % udd = zeros(size(u));
            % udd(2:2:end) = -4;
            % u2Columns = Tri3Element.vector_to_matrix(u); 
            % ud2Columns = Tri3Element.vector_to_matrix(ud);
            % udd2Columns = Tri3Element.vector_to_matrix(udd);
            % uElement = u2Columns(self.nodeIDs,:);
            % udElement = ud2Columns(self.nodeIDs,:);
            % uddElement = udd2Columns(self.nodeIDs,:);
            % nodesPos = self.nodes + uElement; % position = base node position + displacement
            % 
            % % convert back to 6 x 1 (x y x y x y)
            % nodesPos = reshape(nodesPos',[],1);  
            % uElement = reshape(uElement.',[],1);
            % udElement = reshape(udElement',[],1);
            % uddElement = reshape(uddElement',[],1);
            % 
            % T=tensor(T);

            % F = double(ttv(ttv(ttv(T,uddElement,2),nodesPos,2),nodesPos,2) + ttv(ttv(ttv(T,udElement,2),udElement,2),nodesPos,2)+ ttv(ttv(ttv(T,udElement,2),nodesPos,2),udElement,2));
            % test = 2;
            T=double(T);
        end
       
        % HELPER FUNCTIONS ________________________________________________
        function [r11, r12, r21, r22] = normal_vec_rot_matrix(~, startNodePos, endNodePos, nextNodePos)
            n = [0 1; -1 0]*[endNodePos(1)-startNodePos(1); endNodePos(2)-startNodePos(2)];  % 90° rotation clockwise
            n = n/norm(n); % normalizing n
            dotProd = dot(n,[nextNodePos(1)-endNodePos(1); nextNodePos(2)-endNodePos(2)]);
            if dotProd >= 0
                r11 = 0; 
                r12 = -1;
                r21 = 1;
                r22 = 0;
            else
                r11 = 0;
                r12 = 1;
                r21 = -1;
                r22 = 0;
            end
        end 

        function inv = normal_vector_inversion(~, startNodePos, endNodePos, nextNodePos)
            % _____________________________________________________________
            % Returns returns -1 if the normal vector computed with the
            % rotation matrix R=[0 1;-1 0] needs to be inverted Returns 1 
            % if this is not the case.
            % _____________________________________________________________
            n = [0 1; -1 0]*[endNodePos(1)-startNodePos(1); endNodePos(2)-startNodePos(2)];  % 90° rotation clockwise
            n = n/norm(n); % normalizing n
            dotProd = dot(n,[nextNodePos(1)-endNodePos(1); nextNodePos(2)-endNodePos(2)]);

            if dotProd >= 0
                inv = -1; % inversion is needed if n points toward the element
            else
                inv = 1;
            end
        end
        
        function n = normal_vector(~, startNodePos, endNodePos, nextNodePos)
            % _____________________________________________________________
            % Returns a vector of unit length normal to the surface
            % defined by `startNodePos' and enNodePos', pointing
            % outward the considered element
            % _____________________________________________________________
            n = [0 1; -1 0]*[endNodePos(1)-startNodePos(1); endNodePos(2)-startNodePos(2)];  % 90° rotation clockwise
            n = n/norm(n); % normalizing n
            dotProd = dot(n,[nextNodePos(1)-endNodePos(1); nextNodePos(2)-endNodePos(2)]);
            if dotProd >= 0
                n = -n; % inverting the direction of n if it points toward the element
            end
        end

        function [startNode, endNode, nextNode] = get_node_from_face(~, face)
            % _____________________________________________________________
            % Returns the nodes' indexes from a give face. The
            % ordering is done counterclokewise. `nextNode is the index
            % of the node on the next surface.
            % _____________________________________________________________
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

        function tVector = create_tail_t_vector(nodes)
            % _____________________________________________________________
            % Creates a vector tangential to the tail, pointing toward the
            % head of the fish. The "nodes" matrix should contain the tail
            % node in the first row (the two other nodes being in the
            % element but not excatly at the tail).
            % _____________________________________________________________
            % midpoint nodes 2 and 3
            midpoint = 0.5*(nodes(2,:) + nodes(3,:));

            % tail to midpoint normalised vector
            tVector = (midpoint - nodes(1,:))/norm(midpoint - nodes(1,:));
        end

        function pVector = create_tail_p_vector(nodes)
            % _____________________________________________________________
            % Creates a vector perpendicular to the tail, with a x
            % component pointing toward the head of the fish.
            % The "nodes" matrix should contain the tail node in the first 
            % row (the two other nodes being in the element but not excatly
            % at the tail). 
            % _____________________________________________________________            
            % midpoint nodes 2 and 3
            midpoint = 0.5*(nodes(2,:) + nodes(3,:));
        
            % tail to midpoint normalised vector
            tVector = (midpoint - nodes(1,:))/norm(midpoint - nodes(1,:));
        
            % perpendicular vector, with x component pointing toward the head
            if tVector(2)<=0
                pVector = [tVector(2), -tVector(1)];
            else
                pVector = [-tVector(2), tVector(1)];
            end
        end   

        function M = vector_to_matrix(v)
            % _____________________________________________________________
            % Converts a vector with elements corresponding to x and y
            % directions ([x1; y1; x2; y2, x3; y3 ...]) to a matrix
            % separating the x and y elements into 2 columns ([x1 y1; x2
            % y2; x3 y3];
            % _____________________________________________________________
            xComponents = v(1:2:end);
            yComponents = v(2:2:end);
            M = [xComponents, yComponents];
        end
        
    end
        
end % classdef
    

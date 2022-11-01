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

        % Tensors as a function of ud (early-stage version)
        function T = T1e_func_ud(self, specificFace, vwater, rho)
            % _____________________________________________________________
            % Returns the 1st order tensor (i.e., a vector) stemming from 
            % the drag and the thrust force. These forces are acting on the
            % ``specificFace'' of an element. The forces are divided by 2 
            % and applied to the nodes at of the considered face. 
            % _____________________________________________________________
            [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(1));
            length = norm(self.nodes(endNode,:)-self.nodes(startNode,:));
            [r11,r12,r21,r22] = normal_vec_rot_matrix(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
  
            if specificFace(1) == 1
                T = Te_conf_1_drag(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).'...
                    +Te_conf_1_thrust(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).';
            elseif specificFace(1) == 2
                T = Te_conf_2_drag(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).'...
                    +Te_conf_2_thrust(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).';
            else
                T = Te_conf_3_drag(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).'...
                    +Te_conf_3_thrust(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).';
            end
            
            if specificFace(2) ~= 0
                [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(2));
                length = norm(self.nodes(endNode,:)-self.nodes(startNode,:));
                [r11,r12,r21,r22] = normal_vec_rot_matrix(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
                if specificFace(2) == 1
                    T = T + Te_conf_1_drag(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).'...
                        +Te_conf_1_thrust(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).';
                elseif specificFace(2) == 2
                    T = T + Te_conf_2_drag(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).'...
                        +Te_conf_2_thrust(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).';
                else
                    T = T + Te_conf_3_drag(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).'...
                        +Te_conf_3_thrust(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).';
                end
            end
        end

        function T = T2ue(self, specificFace, vwater, rho)
            % _____________________________________________________________
            % Returns a 2nd order tensor (i.e., a matrix) stemming from 
            % the drag and the thrust force. These forces are acting on the
            % ``specificFace'' of an element. The forces are divided by 2 
            % and applied to the nodes at of the considered face. 
            % _____________________________________________________________
            [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(1));
            length = norm(self.nodes(endNode,:)-self.nodes(startNode,:));
            [r11,r12,r21,r22] = normal_vec_rot_matrix(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
           Tue_conf_1_thrust(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).';
            if specificFace(1) == 1
                T = Tue_conf_1_drag(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).';%...
                    %+Tue_conf_1_thrust(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).';
            elseif specificFace(1) == 2
                T = Tue_conf_2_drag(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).'...
                    +Tue_conf_2_thrust(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).';
            else
                T = Tue_conf_3_drag(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).'...
                    +Tue_conf_3_thrust(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).';
            end
            
            if specificFace(2) ~= 0
                [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(2));
                length = norm(self.nodes(endNode,:)-self.nodes(startNode,:));
                [r11,r12,r21,r22] = normal_vec_rot_matrix(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
                if specificFace(2) == 1
                    T = T + Tue_conf_1_drag(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).'...
                        +Tudote_conf_1_thrust(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).';
                elseif specificFace(2) == 2
                    T = T + Tue_conf_2_drag(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).'...
                        +Tue_conf_2_thrust(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).';
                else
                    T = T + Tue_conf_3_drag(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).'...
                        +Tue_conf_3_thrust(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).';
                end
            end
        end

        function T = T2udote(self, specificFace, vwater, rho)
            % _____________________________________________________________
            % Returns a 2nd order tensor (i.e., a matrix) stemming from 
            % the drag and the thrust force. These forces are acting on the
            % ``specificFace'' of an element. The forces are divided by 2 
            % and applied to the nodes at of the considered face. 
            % _____________________________________________________________
            [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(1));
            length = norm(self.nodes(endNode,:)-self.nodes(startNode,:));
            [r11,r12,r21,r22] = normal_vec_rot_matrix(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
  
            if specificFace(1) == 1
                T = Tudote_conf_1_drag(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).'...
                    +Tudote_conf_1_thrust(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).';
            elseif specificFace(1) == 2
                T = Tudote_conf_2_drag(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).'...
                    +Tudote_conf_2_thrust(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).';
            else
                T = Tudote_conf_3_drag(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).'...
                    +Tudote_conf_3_thrust(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).';
            end
            
            if specificFace(2) ~= 0
                [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(2));
                length = norm(self.nodes(endNode,:)-self.nodes(startNode,:));
                [r11,r12,r21,r22] = normal_vec_rot_matrix(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
                if specificFace(2) == 1
                    T = T + Tudote_conf_1_drag(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).'...
                        +Tudote_conf_1_thrust(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).';
                elseif specificFace(2) == 2
                    T = T + Tudote_conf_2_drag(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).'...
                        +Tudote_conf_2_thrust(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).';
                else
                    T = T + Tudote_conf_3_drag(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).'...
                        +Tudote_conf_3_thrust(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).';
                end
            end
        end

        function T = T3uue(self, Ve, specificFace, vwater, rho)
            % _____________________________________________________________
            % Returns a 2nd order tensor (i.e., a matrix) stemming from 
            % the drag and the thrust force. These forces are acting on the
            % ``specificFace'' of an element. The forces are divided by 2 
            % and applied to the nodes at of the considered face. 
            % _____________________________________________________________
            [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(1));
            length = norm(self.nodes(endNode,:)-self.nodes(startNode,:));
            [r11,r12,r21,r22] = normal_vec_rot_matrix(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
     
            if specificFace(1) == 1
                T = Tuue_conf_1_drag(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0)...
                    +Tuue_conf_1_thrust(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0);
            elseif specificFace(1) == 2
                T = Tuue_conf_1_drag(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0)...
                    +Tuue_conf_2_thrust(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0);
            else
                T = Tuue_conf_3_drag(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0)...
                    +Tuue_conf_3_thrust(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0);
            end
            
            if specificFace(2) ~= 0
                [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(2));
                length = norm(self.nodes(endNode,:)-self.nodes(startNode,:));
                [r11,r12,r21,r22] = normal_vec_rot_matrix(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
                if specificFace(2) == 1
                    T = T + Tuue_conf_1_drag(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0)...
                        +Tudote_conf_1_thrust(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0);
                elseif specificFace(2) == 2
                    T = T + Tuue_conf_2_drag(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0)...
                        +Tuue_conf_2_thrust(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0);
                else
                    T = T + Tuue_conf_3_drag(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0)...
                        +Tuue_conf_3_thrust(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0);
                end
            end
            T = tensor(T);
            T = ttm(T,Ve.',1);
            size(T)
        end

        % Tensors in their final form
        function T = Te1(self, specificFace, vwater, rho)
            % _____________________________________________________________
            % Returns the 1st order tensor (i.e., a vector) stemming from 
            % the drag and the thrust force. These forces are acting on the
            % ``specificFace'' of an element. The forces are divided by 2 
            % and applied to the nodes at of the considered face. 
            % _____________________________________________________________
            [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(1));
            length = norm(self.nodes(endNode,:)-self.nodes(startNode,:));
            [r11,r12,r21,r22] = normal_vec_rot_matrix(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
  
            if specificFace(1) == 1
                T = T1_conf1(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2)).';      
            elseif specificFace(1) == 2
                T = T1_conf2(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2)).';   
            else
                T = T1_conf3(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2)).';   
            end
            
            if specificFace(2) ~= 0
                [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(2));
                length = norm(self.nodes(endNode,:)-self.nodes(startNode,:));
                [r11,r12,r21,r22] = normal_vec_rot_matrix(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
                
                if specificFace(2) == 1
                    T = T1_conf1(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2)).';      
                elseif specificFace(1) == 2
                    T = T1_conf2(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2)).';   
                else
                    T = T1_conf3(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2)).';   
                end
            end

        end

        function T = Te2(self, specificFace, vwater, rho)
            % _____________________________________________________________
            % Returns one of the 2nd order tensors (i.e., a matrix) stemming from 
            % the drag and the thrust force. These forces are acting on the
            % ``specificFace'' of an element. The forces are divided by 2 
            % and applied to the nodes at of the considered face. 
            % _____________________________________________________________
            [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(1));
            length = norm(self.nodes(endNode,:)-self.nodes(startNode,:));
            [r11,r12,r21,r22] = normal_vec_rot_matrix(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
  
            if specificFace(1) == 1
                T = T2_conf1(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));      
            elseif specificFace(1) == 2
                T = T2_conf2(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
            else
                T = T2_conf3(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
            end
            
            if specificFace(2) ~= 0
                [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(2));
                length = norm(self.nodes(endNode,:)-self.nodes(startNode,:));
                [r11,r12,r21,r22] = normal_vec_rot_matrix(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
                
                if specificFace(2) == 1
                    T = T2_conf1(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));      
                elseif specificFace(1) == 2
                    T = T2_conf2(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
                else
                    T = T2_conf3(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
                end
            end

        end

        
        function T = Teu2(self, specificFace, vwater, rho)
            % _____________________________________________________________
            % Returns one of the 2nd order tensors (i.e., a matrix) stemming from 
            % the drag and the thrust force. These forces are acting on the
            % ``specificFace'' of an element. The forces are divided by 2 
            % and applied to the nodes at of the considered face. 
            % _____________________________________________________________
            [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(1));
            length = norm(self.nodes(endNode,:)-self.nodes(startNode,:));
            [r11,r12,r21,r22] = normal_vec_rot_matrix(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
  
            if specificFace(1) == 1
                T = Tu2_conf1(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));      
            elseif specificFace(1) == 2
                T = Tu2_conf2(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
            else
                T = Tu2_conf3(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
            end
            
            if specificFace(2) ~= 0
                [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(2));
                length = norm(self.nodes(endNode,:)-self.nodes(startNode,:));
                [r11,r12,r21,r22] = normal_vec_rot_matrix(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
                
                if specificFace(2) == 1
                    T = Tu2_conf1(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));      
                elseif specificFace(1) == 2
                    T = Tu2_conf2(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
                else
                    T = Tu2_conf3(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
                end
            end

        end

        function T = Teu3(self, specificFace, vwater, rho)
            % _____________________________________________________________
            % Returns one of the 3rd order tensors stemming from 
            % the drag and the thrust force. These forces are acting on the
            % ``specificFace'' of an element. The forces are divided by 2 
            % and applied to the nodes at of the considered face. 
            % _____________________________________________________________
            [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(1));
            length = norm(self.nodes(endNode,:)-self.nodes(startNode,:));
            [r11,r12,r21,r22] = normal_vec_rot_matrix(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
  
            if specificFace(1) == 1
                T = Tu3_conf1(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));      
            elseif specificFace(1) == 2
                T = Tu3_conf2(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
            else
                T = Tu3_conf3(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
            end
            
            if specificFace(2) ~= 0
                [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(2));
                length = norm(self.nodes(endNode,:)-self.nodes(startNode,:));
                [r11,r12,r21,r22] = normal_vec_rot_matrix(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
                
                if specificFace(2) == 1
                    T = Tu3_conf1(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));      
                elseif specificFace(1) == 2
                    T = Tu3_conf2(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
                else
                    T = Tu3_conf3(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
                end
            end

            T = permute(T,[3 1 2]);

        end

        function T = Teudot2(self, specificFace, vwater, rho)
            % _____________________________________________________________
            % Returns one of the 2nd order tensors (i.e., a matrix) stemming from 
            % the drag and the thrust force. These forces are acting on the
            % ``specificFace'' of an element. The forces are divided by 2 
            % and applied to the nodes at of the considered face. 
            % _____________________________________________________________
            [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(1));
            length = norm(self.nodes(endNode,:)-self.nodes(startNode,:));
            [r11,r12,r21,r22] = normal_vec_rot_matrix(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
  
            if specificFace(1) == 1
                T = Tudot2_conf1(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));      
            elseif specificFace(1) == 2
                T = Tudot2_conf2(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
            else
                T = Tudot2_conf3(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
            end
            
            if specificFace(2) ~= 0
                [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(2));
                length = norm(self.nodes(endNode,:)-self.nodes(startNode,:));
                [r11,r12,r21,r22] = normal_vec_rot_matrix(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
                
                if specificFace(2) == 1
                    T = Tudot2_conf1(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));      
                elseif specificFace(1) == 2
                    T = Tudot2_conf2(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
                else
                    T = Tudot2_conf3(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
                end
            end

        end

        function T = Teudot3(self, specificFace, vwater, rho)
            % _____________________________________________________________
            % Returns one of the 2nd order tensors (i.e., a matrix) stemming from 
            % the drag and the thrust force. These forces are acting on the
            % ``specificFace'' of an element. The forces are divided by 2 
            % and applied to the nodes at of the considered face. 
            % _____________________________________________________________
            [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(1));
            length = norm(self.nodes(endNode,:)-self.nodes(startNode,:));
            [r11,r12,r21,r22] = normal_vec_rot_matrix(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
  
            if specificFace(1) == 1
                T = Tudot3_conf1(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));      
            elseif specificFace(1) == 2
                T = Tudot3_conf2(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
            else
                T = Tudot3_conf3(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
            end
            
            if specificFace(2) ~= 0
                [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(2));
                length = norm(self.nodes(endNode,:)-self.nodes(startNode,:));
                [r11,r12,r21,r22] = normal_vec_rot_matrix(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
                
                if specificFace(2) == 1
                    T = Tudot3_conf1(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));      
                elseif specificFace(1) == 2
                    T = Tudot3_conf2(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
                else
                    T = Tudot3_conf3(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
                end
            end

        end
        
        
        function T = Teuu3(self, specificFace, vwater, rho)
            % _____________________________________________________________
            % Returns one of the 3rd order tensors stemming from 
            % the drag and the thrust force. These forces are acting on the
            % ``specificFace'' of an element. The forces are divided by 2 
            % and applied to the nodes at of the considered face. 
            % _____________________________________________________________
            [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(1));
            length = norm(self.nodes(endNode,:)-self.nodes(startNode,:));
            [r11,r12,r21,r22] = normal_vec_rot_matrix(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
            if specificFace(1) == 1
                T = Tuu3_conf1(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));      
            elseif specificFace(1) == 2
                T = Tuu3_conf2(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
            else
                T = Tuu3_conf3(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
            end
            
            if specificFace(2) ~= 0
                [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(2));
                length = norm(self.nodes(endNode,:)-self.nodes(startNode,:));
                [r11,r12,r21,r22] = normal_vec_rot_matrix(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
                
                if specificFace(2) == 1
                    T = Tuu3_conf1(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));      
                elseif specificFace(1) == 2
                    T = Tuu3_conf2(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
                else
                    T = Tuu3_conf3(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
                end
            end

            T = permute(T,[3 1 2]);
        end

        function T = Teuudot3(self, specificFace, vwater, rho)
            % _____________________________________________________________
            % Returns one of the 3rd order tensors stemming from 
            % the drag and the thrust force. These forces are acting on the
            % ``specificFace'' of an element. The forces are divided by 2 
            % and applied to the nodes at of the considered face. 
            % _____________________________________________________________
            [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(1));
            length = norm(self.nodes(endNode,:)-self.nodes(startNode,:));
            [r11,r12,r21,r22] = normal_vec_rot_matrix(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
            if specificFace(1) == 1
                T = Tuudot3_conf1(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));      
            elseif specificFace(1) == 2
                T = Tuudot3_conf2(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
            else
                T = Tuudot3_conf3(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
            end
            
            if specificFace(2) ~= 0
                [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(2));
                length = norm(self.nodes(endNode,:)-self.nodes(startNode,:));
                [r11,r12,r21,r22] = normal_vec_rot_matrix(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
                
                if specificFace(2) == 1
                    T = Tuudot3_conf1(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));      
                elseif specificFace(1) == 2
                    T = Tuudot3_conf2(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
                else
                    T = Tuudot3_conf3(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
                end
            end

            T = permute(T,[3 1 2]);
        end

        function T = Teudotudot3(self, specificFace, vwater, rho)
            % _____________________________________________________________
            % Returns one of the 3rd order tensors stemming from 
            % the drag and the thrust force. These forces are acting on the
            % ``specificFace'' of an element. The forces are divided by 2 
            % and applied to the nodes at of the considered face. 
            % _____________________________________________________________
            [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(1));
            length = norm(self.nodes(endNode,:)-self.nodes(startNode,:));
            [r11,r12,r21,r22] = normal_vec_rot_matrix(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
            if specificFace(1) == 1
                T = Tudotudot3_conf1(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));      
            elseif specificFace(1) == 2
                T = Tudotudot3_conf2(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
            else
                T = Tudotudot3_conf3(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
            end
            
            if specificFace(2) ~= 0
                [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(2));
                length = norm(self.nodes(endNode,:)-self.nodes(startNode,:));
                [r11,r12,r21,r22] = normal_vec_rot_matrix(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
                
                if specificFace(2) == 1
                    T = Tudotudot3_conf1(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));      
                elseif specificFace(1) == 2
                    T = Tudotudot3_conf2(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
                else
                    T = Tudotudot3_conf3(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2));   
                end
            end

            T = permute(T,[3 1 2]);
          end

        % Forces for early-stage tests
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

        function F = thrust_force(self, specificFace, vwater, rho)
            % _____________________________________________________________
            % Returns the thrust force F_drag acting the face
            % ``specificFace'' of an element. The force is divided by 2 and
            % applied to the nodes at of the considered face. Use of 1st
            % order tensors.
            % _____________________________________________________________
            [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(1));
            length = norm(self.nodes(endNode,:)-self.nodes(startNode,:));
            [r11,r12,r21,r22] = normal_vec_rot_matrix(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
           
            if specificFace(1) == 1
                F = Te_conf_1_thrust(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).';
            elseif specificFace(1) == 2
                F = Te_conf_2_thrust(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).';
            else
                F = Te_conf_3_thrust(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).';
            end
             %disp(F)
            if specificFace(2) ~= 0
                [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(2));
                length = norm(self.nodes(endNode,:)-self.nodes(startNode,:));
                [r11,r12,r21,r22] = normal_vec_rot_matrix(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
                if specificFace(2) == 1
                    F = F + Te_conf_1_thrust(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).';
                elseif specificFace(2) == 2
                    F = F + Te_conf_2_thrust(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).';
                else
                    F = F + Te_conf_3_thrust(length,rho,vwater(1),vwater(2),r11,r12,r21,r22,self.nodes(1,1),self.nodes(1,2),self.nodes(2,1),self.nodes(2,2),self.nodes(3,1),self.nodes(3,2),0,0,0,0,0,0).';
                end
            end
            
            
        end

        function F = drag_force_full(self, specificFace, vwater, rho)
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
                n = normal_vector(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));  
                length = norm(self.nodes(endNode,:)-self.nodes(startNode,:));
                A = length; % in 2D, area proportional force is length proportional
                ns = 2; % 2D surfaces of TRI3 are composed of 2 elements
                vrel = vwater - 1/ns*[0;0]; %(self.nodes(startNode,:) + self.nodes(endNode,:)); % need to be replaced by velocity
                d = vrel/norm(vrel);
                Cd = 2*(n.'*d)^2;
                
                % Drag force
                Fdrag = 1/2*rho*A*Cd*norm(vrel)^2*d;
            end 
            
        end

        function F = thrust_force_full(self, specificFace, vwater, rho)
            % _____________________________________________________________
            % Returns the thrust force F_drag acting the face
            % ``specificFace'' of an element. The force is divided by 2 and
            % applied to the nodes at of the considered face.
            % _____________________________________________________________
            F = sparse(self.nelDOFs,1);
            [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(1));

            % Computes thrust force
            Fthrust = compute_thrust_force(startNode, endNode, nextNode);
            % Applies thrust force on nodes
            F = apply_force(self, F, Fthrust, startNode, endNode);

            if specificFace(2) ~= 0
                [startNode, endNode, nextNode] = get_node_from_face(self,specificFace(2));

                % Computes thrust force
                Fthrust = compute_thrust_force(startNode, endNode, nextNode);
                % Apply thrust force on nodes
                F = apply_force(self, F, Fthrust, startNode, endNode);

            end

            function Fthrust = compute_thrust_force(startNode, endNode, nextNode)
                % _________________________________________________________
                % Computes the thrust force vector for a specific surface
                %__________________________________________________________
                % Terms needed for the thrust force - REPLACE VREL (to do) !!
                n = normal_vector(self, self.nodes(startNode,:), self.nodes(endNode,:), self.nodes(nextNode,:));
                length = norm(self.nodes(endNode,:)-self.nodes(startNode,:));
                A = length; % in 2D, area proportional force is length proportional
                ns = 2; % 2D surfaces of TRI3 are composed of 2 elements
                vrel = vwater - 1/ns*[0;0]; %(self.nodes(startNode,:) + self.nodes(endNode,:)); % need to be replaced by velocity
                d = vrel/norm(vrel);
                Ct = 1/3*((acos(n.'*d))^2-pi^2/4);
                
                % Thrust force
                Fthrust = -1/2*A*rho*Ct*norm(vrel)^2*n;
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
        
        % Useful helper functions
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
        
    end
        
end % classdef
    

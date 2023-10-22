classdef ReducedAssembly < Assembly
    properties
        % Already contains properties in the Assembly class
        
        V       % Reduction basis
    end

    methods
        function self = ReducedAssembly(Mesh,V)
            % Assembly (superclass) constructor
            self@Assembly(Mesh);
            
            % set reduction basis
            self.V = V;
        end

        function set.V(self,V)
            if size(V,1) ~= self.Mesh.nDOFs %#ok<*MCSUP>
                warning(['Reduction basis size incorrect: should have ' ...
                    num2str(self.Mesh.nDOFs) ' rows'])
            end
            self.V = V;
        end

        function [f] = uniform_body_force(self,varargin)
            m = size(self.V,2);
            f = zeros(m,1);
            V = self.V;             %#ok<*PROP>
            [elementWeights,~] = self.parse_inputs(varargin{:});
            
            % extracting elements with nonzero weights
            elementSet = find(elementWeights);
            
            % extracting domain elements from the supplied set
            domainElements = ~[self.Mesh.Elements(elementSet).isBoundary];
            
            Elements = self.Mesh.Elements;
            parfor j = elementSet(domainElements)
                thisElement = Elements(j).Object;
                index = thisElement.iDOFs;          
                Ve = V(index,:);
                fe = thisElement.uniformBodyForce;
                f = f + elementWeights(j) * (Ve.' * fe);
            end
        end


        function [K] = matrix(self,elementMethodName,varargin)
            % This function assembles a generic finite element matrix from
            % its element level counterpart.
            % elementMethodName is a string input containing the name of
            % the method that returns the element level matrix Ke.
            % For this to work, a method named elementMethodName which
            % returns the appropriate matrix must be defined for all
            % element types in the FE Mesh.            
            % NOTE: it is assumed that the input arguments are provided in
            % the full (unreduced) system
            
            m = size(self.V,2);
            K = zeros(m,m);
            Elements = self.Mesh.Elements;
            V = self.V;
            % parsing element weights
            [elementWeights,inputs] = self.parse_inputs(varargin{:});
            
            % extracting elements with nonzero weights
            elementSet = find(elementWeights);
            
            % Computing element level contributions
            parfor j = elementSet
                thisElement = Elements(j).Object;
                index = thisElement.iDOFs;          
                Ve = V(index,:);
                Ke = thisElement.(elementMethodName)(inputs{:});
                K = K + elementWeights(j) * (Ve.' * Ke * Ve);
            end
        end


        function [K] = matrix_actuation(self,elementMethodName,varargin)
            % This function assembles a generic finite element matrix from
            % its element level counterpart. It uses a for loop instead of
            % a parfor loop to allow the elementSet to be a subpart of the
            % whole elements set.
            % elementMethodName is a string input containing the name of
            % the method that returns the element level matrix Ke.
            % For this to work, a method named elementMethodName which
            % returns the appropriate matrix must be defined for all
            % element types in the FE Mesh.            
            % NOTE: it is assumed that the input arguments are provided in
            % the full (unreduced) system
            
            m = size(self.V,2);
            K = zeros(m,m);
            Elements = self.Mesh.Elements;
            V = self.V;
            % parsing element weights
            [elementWeights,inputs] = self.parse_inputs(varargin{:});
            
            % extracting elements with nonzero weights
            elementSet = find(elementWeights);
            
            % Computing element level contributions
            for j = elementSet
                thisElement = Elements(j).Object;
                index = thisElement.iDOFs;          
                Ve = V(index,:);
                Ke = thisElement.(elementMethodName)(inputs{:});
                K = K + elementWeights(j) * (Ve.' * Ke * Ve);
            end
        end

        function [K] = matrix_actuation_PROM(self,elementMethodName,U,varargin)
            % This function assembles a generic finite element matrix from
            % its element level counterpart. It uses a for loop instead of
            % a parfor loop to allow the elementSet to be a subpart of the
            % whole elements set.
            % elementMethodName is a string input containing the name of
            % the method that returns the element level matrix Ke.
            % For this to work, a method named elementMethodName which
            % returns the appropriate matrix must be defined for all
            % element types in the FE Mesh.            
            % NOTE: it is assumed that the input arguments are provided in
            % the full (unreduced) system
            
            m = size(self.V,2);
            md = size(U,2);
            K = zeros(m,md);
            Elements = self.Mesh.Elements;
            V = self.V;
            % parsing element weights
            [elementWeights,inputs] = self.parse_inputs(varargin{:});
            
            % extracting elements with nonzero weights
            elementSet = find(elementWeights);
            
            % Computing element level contributions
            for j = elementSet
                thisElement = Elements(j).Object;
                index = thisElement.iDOFs;          
                Ve = V(index,:);
                Ue = U(index,:);
                Ke = thisElement.(elementMethodName)(inputs{:});
                K = K + elementWeights(j) * (Ve.' * Ke * Ue);
            end
        end

        function [K] = matrix_skin(self,elementMethodName,varargin)
            % This function assembles a finite element matrix from
            % its element level counterpart. The method allows to pass
            % extra argument to access and work with the skin
            % elements/nodes of the structure.
            % elementMethodName is a string input containing the name of
            % the method that returns the element level matrix Ke.
            % For this to work, a method named elementMethodName which
            % returns the appropriate matrix must be defined for all
            % element types in the FE Mesh.            
            % NOTE: it is assumed that the input arguments are provided in
            % the full (unreduced) system
            
            m = size(self.V,2);
            K = zeros(m,m);
            Elements = self.Mesh.Elements;
            V = self.V;
            % parsing element weights
            [elementWeights,inputs] = self.parse_inputs(varargin{:});
            
            % extracting elements with nonzero weights
            elementSet = find(elementWeights);
            
            % Computing element level contributions
            for j = elementSet
                thisElement = Elements(j).Object;
                index = thisElement.iDOFs;          
                Ve = V(index,:);
                Ke = thisElement.(elementMethodName)(inputs{1}(j,:),inputs{2});
                K = K + elementWeights(j) * (Ve.' * Ke * Ve);
            end
        end

        function [K] = matrix_skin_PROM(self,elementMethodName,U,varargin)
            % This function assembles a finite element matrix from
            % its element level counterpart. The method allows to pass
            % extra argument to access and work with the skin
            % elements/nodes of the structure.
            % elementMethodName is a string input containing the name of
            % the method that returns the element level matrix Ke.
            % For this to work, a method named elementMethodName which
            % returns the appropriate matrix must be defined for all
            % element types in the FE Mesh.            
            % NOTE: it is assumed that the input arguments are provided in
            % the full (unreduced) system
            
            m = size(self.V,2);
            md = size(U,2);
            K = zeros(m,md);
            Elements = self.Mesh.Elements;
            V = self.V;
            % parsing element weights
            [elementWeights,inputs] = self.parse_inputs(varargin{:});
            
            % extracting elements with nonzero weights
            elementSet = find(elementWeights);
            
            % Computing element level contributions
            for j = elementSet
                thisElement = Elements(j).Object;
                index = thisElement.iDOFs;          
                Ve = V(index,:);
                Ue = U(index,:);
                Ke = thisElement.(elementMethodName)(inputs{1}(j,:),inputs{2}, inputs{3}, inputs{4});
                K = K + elementWeights(j) * (Ve.' * Ke * Ue);
            end
        end
        

        function [f] = vector(self,elementMethodName,varargin)
            % This function assembles a generic finite element vector from
            % its element level counterpart.
            % elementMethodName is a string input containing the name of
            % the method that returns the element level vector Fe.
            % For this to work, a method named elementMethodName which
            % returns the appropriate vector must be defined for all
            % element types in the FE Mesh.            
            % NOTE: it is assumed that the input arguments are provided in
            % the full (unreduced) system

            m = size(self.V,2);
            f = zeros(m,1);

            Elements = self.Mesh.Elements;
            V = self.V;
            % parsing element weights
            [elementWeights,inputs] = self.parse_inputs(varargin{:});
            
            % extracting elements with nonzero weights
            elementSet = find(elementWeights);
            
            % Computing element level contributions
            parfor j = elementSet 
                thisElement = Elements(j).Object;
                index = thisElement.iDOFs;          
                Ve = V(index,:);
                fe = thisElement.(elementMethodName)(inputs{:});
                f = f + elementWeights(j) * (Ve.' * fe);
            end
        end

        function [f] = vector_actuation(self,elementMethodName,varargin)
            % This function assembles a generic finite element vector from
            % its element level counterpart. It uses a for loop instead of
            % a parfor loop to allow the elementSet to be a subpart of the
            % whole elements set.
            % elementMethodName is a string input containing the name of
            % the method that returns the element level vector Fe.
            % For this to work, a method named elementMethodName which
            % returns the appropriate vector must be defined for all
            % element types in the FE Mesh.            
            % NOTE: it is assumed that the input arguments are provided in
            % the full (unreduced) system

            m = size(self.V,2);
            f = zeros(m,1);

            Elements = self.Mesh.Elements;
            V = self.V;
            % parsing element weights
            [elementWeights,inputs] = self.parse_inputs(varargin{:});
            
            % extracting elements with nonzero weights
            elementSet = find(elementWeights);
            
            % Computing element level contributions
            for j = elementSet 
                thisElement = Elements(j).Object;
                index = thisElement.iDOFs;          
                Ve = V(index,:);
                fe = thisElement.(elementMethodName)(inputs{:});
                f = f + elementWeights(j) * (Ve.' * fe);
            end
        end

        function [f] = vector_skin(self,elementMethodName,varargin)
            % This function assembles a finite element vector from
            % its element level counterpart. The method allows to pass
            % extra argument to access and work with the skin
            % elements/nodes of the structure.
            % elementMethodName is a string input containing the name of
            % the method that returns the element level vector Fe.
            % For this to work, a method named elementMethodName which
            % returns the appropriate vector must be defined for all
            % element types in the FE Mesh.            
            % NOTE: it is assumed that the input arguments are provided in
            % the full (unreduced) system

            m = size(self.V,2);
            f = zeros(m,1);

            Elements = self.Mesh.Elements;
            V = self.V;
            % parsing element weights
            [elementWeights,inputs] = self.parse_inputs(varargin{:});
            
            % extracting elements with nonzero weights
            elementSet = find(elementWeights);
            
            % Computing element level contributions
            for j = elementSet %was a parfor
                thisElement = Elements(j).Object;
                index = thisElement.iDOFs;          
                Ve = V(index,:);
                fe = thisElement.(elementMethodName)(inputs{1}(j,:),inputs{2});
                f = f + elementWeights(j) * (Ve.' * fe);
            end
        end
        
        
        function [K, f] = matrix_and_vector(self,elementMethodName, varargin)
            % This function assembles a generic finite element matrix and
            % vector from its element level counterpart.
            % elementMethodName is a string input containing the name of
            % the method that returns the element level matrix Ke.
            % For this to work, a method named elementMethodName which
            % returns the appropriate matrix must be defined for all
            % element types in the FE Mesh.            
            % NOTE: it is assumed that the input arguments are provided in
            % the full (unreduced) system
            
            m = size(self.V,2);
            
            K = zeros(m,m);
            f = zeros(m,1);
            Elements = self.Mesh.Elements;
            V = self.V;
            % parsing element weights
            [elementWeights,inputs] = self.parse_inputs(varargin{:});
            
            % extracting elements with nonzero weights
            elementSet = find(elementWeights);
            
            % Computing element level contributions
            parfor j = elementSet
                thisElement = Elements(j).Object;
                index = thisElement.iDOFs;          
                Ve = V(index,:);
                [Ke, fe] = thisElement.(elementMethodName)(inputs{:});
                K = K + elementWeights(j) * (Ve.' * Ke * Ve);
                f = f + elementWeights(j) * (Ve.' * fe);
            end
        end
        
        
        function [T] = tensor(self,elementMethodName,SIZE,sumDIMS,mode,varargin)
            % This function assembles a generic finite element vector from
            % its element level counterpart.
            % elementMethodName is a string input containing the name of
            % the method that returns the element level vector Fe.
            % For this to work, a method named elementMethodName which
            % returns the appropriate vector must be defined for all
            % element types in the FE Mesh.            
            % NOTE: in this function, we reduce all the dimensions with 
            % the same reduction basis
            
            m = size(self.V,2);
            [~,I] = find(SIZE == m);
            T = tenzeros(SIZE);
            Elements = self.Mesh.Elements;
            V = self.V;                %#ok<*PROPLC>
            
            % parsing element weights
            [elementWeights,inputs] = self.parse_inputs(varargin{:});
            
            % extracting elements with nonzero weights
            elementSet = find(elementWeights);
            
            % Computing element level contributions

            for j = elementSet
                thisElement = Elements(j).Object;
                index = thisElement.iDOFs;          
                Ve = V(index,:); %#ok<*PFBNS>
                if strcmpi(mode,'ELP') % toggle Element-Level projection
                    [Te, ~] = thisElement.(elementMethodName)([Ve,inputs{:}]);
                    T = T + Te;
                else
                    [Te, ~] = thisElement.(elementMethodName)(inputs{:});
                    % transform tensor
                    Vcell = cell(ndims(Te),1);
                    Vcell(:) = {Ve.'};
                    T = T + elementWeights(j) * ttm(Te, Vcell,I);
                end
            end
            
            [subs, T] = sparsify(T,[],sumDIMS);
            T = sptensor(subs, T, SIZE);
        end

        
        function [T] = tensor_spine_momentum_xx(self,elementMethodName,SIZE,varargin)
            % 
            
            T = tenzeros(SIZE);
            Elements = self.Mesh.Elements;
            V = self.V;              

            % parsing element weights
            [elementWeights,inputs] = self.parse_inputs(varargin{:});
            
            % extracting elements with nonzero weights
            elementSet = find(elementWeights);
            
            % Computing element level contributions

            for j = elementSet
                thisElement = Elements(j).Object;
                index = thisElement.iDOFs;   
                x0 = reshape(thisElement.nodes.',[],1);
                Ve = V(index,:); %#ok<*PFBNS>
                
                Te = thisElement.(elementMethodName)(inputs{1}(j,:),inputs{2}(j), inputs{3});
                Ter = einsum('iI,ijkl,jJ,kK,lL->IJKL',Ve,Te,Ve,x0,x0);  % will be of size m x m, as the two last dimensions of 1 will not appear in Ter
 
                T = T + Ter;

            end
            T = tensor(T);
            
        end

        function [T] = tensor_spine_momentum_xV(self,elementMethodName,SIZE,varargin)
            % 
            
            T = tenzeros(SIZE);
            Elements = self.Mesh.Elements;
            V = self.V; 
            m = size(V,2);

            % parsing element weights
            [elementWeights,inputs] = self.parse_inputs(varargin{:});
            
            % extracting elements with nonzero weights
            elementSet = find(elementWeights);
            
            % Computing element level contributions

            for j = elementSet
                thisElement = Elements(j).Object;
                index = thisElement.iDOFs;
                x0 = reshape(thisElement.nodes.',[],1);
                Ve = V(index,:); %#ok<*PFBNS>
                
                Te = thisElement.(elementMethodName)(inputs{1}(j,:),inputs{2}(j), inputs{3});
                Ter = einsum('iI,ijkl,jJ,kK,lL->IJKL',Ve,Te,Ve,x0,Ve); % will be of size m x m x 1 x m
                
                if m ~= 0
                    Ter = ttv(tensor(Ter),1,3);   % bring it in the form m x m x m
                end
                
                T = T + Ter;

            end
            T = tensor(T);
            
        end

        function [T] = tensor_spine_momentum_Vx(self,elementMethodName,SIZE,varargin)
            % 
            
            T = tenzeros(SIZE);
            Elements = self.Mesh.Elements;
            V = self.V;              

            % parsing element weights
            [elementWeights,inputs] = self.parse_inputs(varargin{:});
            
            % extracting elements with nonzero weights
            elementSet = find(elementWeights);
            
            % Computing element level contributions

            for j = elementSet
                thisElement = Elements(j).Object;
                index = thisElement.iDOFs;   
                x0 = reshape(thisElement.nodes.',[],1);
                Ve = V(index,:); %#ok<*PFBNS>
                
                Te = thisElement.(elementMethodName)(inputs{1}(j,:),inputs{2}(j), inputs{3});
                Ter = einsum('iI,ijkl,jJ,kK,lL->IJKL',Ve,Te,Ve,Ve,x0); % will be of size m x m x m, as the forth dimension of 1 will not appear in Ter

                T = T + Ter;

            end
            T = tensor(T);
            
        end

        function [T] = tensor_spine_momentum_VV(self,elementMethodName,SIZE,varargin)
            % 
            
            T = tenzeros(SIZE);
            Elements = self.Mesh.Elements;
            V = self.V;              

            % parsing element weights
            [elementWeights,inputs] = self.parse_inputs(varargin{:});
            
            % extracting elements with nonzero weights
            elementSet = find(elementWeights);
            
            % Computing element level contributions

            for j = elementSet
                thisElement = Elements(j).Object;
                index = thisElement.iDOFs;          
                Ve = V(index,:); %#ok<*PFBNS>
                
                Te = thisElement.(elementMethodName)(inputs{1}(j,:),inputs{2}(j), inputs{3});
                Ter = einsum('iI,ijkl,jJ,kK,lL->IJKL',Ve,Te,Ve,Ve,Ve);

                T = T + Ter;

            end
            T = tensor(T);
            
        end

        function [T] = tensor_spine_momentum_xx_TET4(self,elementMethodName,SIZE,varargin)
            % 
            
            T = tenzeros(SIZE);
            Elements = self.Mesh.Elements;
            V = self.V;              

            % parsing element weights
            [elementWeights,inputs] = self.parse_inputs(varargin{:});
            
            % extracting elements with nonzero weights
            elementSet = find(elementWeights);
            
            % Computing element level contributions

            for j = elementSet
                thisElement = Elements(j).Object;
                index = thisElement.iDOFs;   
                x0 = reshape(thisElement.nodes.',[],1);
                Ve = V(index,:); %#ok<*PFBNS>
                
                Te = thisElement.(elementMethodName)(inputs{1}(j,:),inputs{2}(j), inputs{3}(j));
                Ter = einsum('iI,ijkl,jJ,kK,lL->IJKL',Ve,Te,Ve,x0,x0);  % will be of size m x m, as the two last dimensions of 1 will not appear in Ter
 
                T = T + Ter;

            end
            T = tensor(T);
            
        end

        function [T] = tensor_spine_momentum_xV_TET4(self,elementMethodName,SIZE,varargin)
            % 
            
            T = tenzeros(SIZE);
            Elements = self.Mesh.Elements;
            V = self.V; 
            m = size(V,2);

            % parsing element weights
            [elementWeights,inputs] = self.parse_inputs(varargin{:});
            
            % extracting elements with nonzero weights
            elementSet = find(elementWeights);
            
            % Computing element level contributions

            for j = elementSet
                thisElement = Elements(j).Object;
                index = thisElement.iDOFs;
                x0 = reshape(thisElement.nodes.',[],1);
                Ve = V(index,:); %#ok<*PFBNS>
                
                Te = thisElement.(elementMethodName)(inputs{1}(j,:),inputs{2}(j), inputs{3}(j));
                Ter = einsum('iI,ijkl,jJ,kK,lL->IJKL',Ve,Te,Ve,x0,Ve); % will be of size m x m x 1 x m
                
                if m ~= 0
                    Ter = ttv(tensor(Ter),1,3);   % bring it in the form m x m x m
                end
                
                T = T + Ter;

            end
            T = tensor(T);
            
        end

        function [T] = tensor_spine_momentum_Vx_TET4(self,elementMethodName,SIZE,varargin)
            % 
            
            T = tenzeros(SIZE);
            Elements = self.Mesh.Elements;
            V = self.V;              

            % parsing element weights
            [elementWeights,inputs] = self.parse_inputs(varargin{:});
            
            % extracting elements with nonzero weights
            elementSet = find(elementWeights);
            
            % Computing element level contributions

            for j = elementSet
                thisElement = Elements(j).Object;
                index = thisElement.iDOFs;   
                x0 = reshape(thisElement.nodes.',[],1);
                Ve = V(index,:); %#ok<*PFBNS>
                
                Te = thisElement.(elementMethodName)(inputs{1}(j,:),inputs{2}(j), inputs{3}(j));
                Ter = einsum('iI,ijkl,jJ,kK,lL->IJKL',Ve,Te,Ve,Ve,x0); % will be of size m x m x m, as the forth dimension of 1 will not appear in Ter

                T = T + Ter;

            end
            T = tensor(T);
            
        end

        function [T] = tensor_spine_momentum_VV_TET4(self,elementMethodName,SIZE,varargin)
            % 
            
            T = tenzeros(SIZE);
            Elements = self.Mesh.Elements;
            V = self.V;              

            % parsing element weights
            [elementWeights,inputs] = self.parse_inputs(varargin{:});
            
            % extracting elements with nonzero weights
            elementSet = find(elementWeights);
            
            % Computing element level contributions

            for j = elementSet
                thisElement = Elements(j).Object;
                index = thisElement.iDOFs;          
                Ve = V(index,:); %#ok<*PFBNS>
                
                Te = thisElement.(elementMethodName)(inputs{1}(j,:),inputs{2}(j), inputs{3}(j));
                Ter = einsum('iI,ijkl,jJ,kK,lL->IJKL',Ve,Te,Ve,Ve,Ve);

                T = T + Ter;

            end
            T = tensor(T);
            
        end

        function [T] = tensor_spine_momentum_TET4(self,elementMethodName,SIZE,varargin)
            % This function assembles a generic finite element vector from
            % its element level counterpart.
            % elementMethodName is a string input containing the name of
            % the method that returns the element level vector Fe.
            % For this to work, a method named elementMethodName which
            % returns the appropriate vector must be defined for all
            % element types in the FE Mesh.            
            % NOTE: in this function, we reduce all the dimensions with 
            % the same reduction basis
            
            T = tenzeros(SIZE);
            Elements = self.Mesh.Elements;
            V = self.V;              

            % parsing element weights
            [elementWeights,inputs] = self.parse_inputs(varargin{:});
            
            % extracting elements with nonzero weights
            elementSet = find(elementWeights);
            
            % Computing element level contributions

            for j = elementSet
                thisElement = Elements(j).Object;
                index = thisElement.iDOFs;          
                Ve = V(index,:); %#ok<*PFBNS>
                
                Te = thisElement.(elementMethodName)(inputs{1}(j,:),inputs{2}(j), inputs{3}(j));
                Ter = einsum('iI,ijkl,jJ,kK,lL->IJKL',Ve,Te,Ve,Ve,Ve);
                % test
                % =einsum('iI,ijkl,jJ,kK,lL->IJKL',Ve',Ter,Ve',Ve',Ve') %
                % map back
                T = T + Ter;

            end
            T = tensor(T);
            
        end
        
        function [T] = tensor_skin(self,elementMethodName,SIZE,sumDIMS,mode,varargin)
            % This function assembles a generic finite element vector from
            % its element level counterpart.
            % elementMethodName is a string input containing the name of
            % the method that returns the element level vector Fe.
            % For this to work, a method named elementMethodName which
            % returns the appropriate vector must be defined for all
            % element types in the FE Mesh.            
            % NOTE: in this function, we reduce all the dimensions with 
            % the same reduction basis
            
            m = size(self.V,2);
            [~,I] = find(SIZE == m);
            T = tenzeros(SIZE);
            Elements = self.Mesh.Elements;
            V = self.V;                %#ok<*PROPLC>
            
            % parsing element weights
            [elementWeights,inputs] = self.parse_inputs(varargin{:});
            
            % extracting elements with nonzero weights
            elementSet = find(elementWeights);
            
            % Computing element level contributions

            for j = elementSet
                thisElement = Elements(j).Object;
                index = thisElement.iDOFs;          
                Ve = V(index,:); %#ok<*PFBNS>
                if strcmpi(mode,'ELP') % toggle Element-Level projection
                    Te = thisElement.(elementMethodName)(inputs{1}(j,:),inputs{2});
                    Ter = einsum('iI,ijk,jJ,kK->IJK',Ve,Te,Ve,Ve);
                    T = T + Ter;
                else
                      msg = 'Only the ELP mode is currently supported for hydrodynamic forces';
                      error(msg)
                end
            end
            T = tensor(T);
            % Next possible improvement: store the tensor as a sparse
            % tensor
            %[subs, T] = sparsify(T,[],sumDIMS);
            %T = sptensor(subs, T, SIZE);
        end

        function [T] = tensor_skin_PROM(self,elementMethodName,U,SIZE,sumDIMS,mode,varargin)
            % This function assembles a generic finite element vector from
            % its element level counterpart.
            % elementMethodName is a string input containing the name of
            % the method that returns the element level vector Fe.
            % For this to work, a method named elementMethodName which
            % returns the appropriate vector must be defined for all
            % element types in the FE Mesh.            
            % NOTE: in this function, we reduce all the dimensions with 
            % the same reduction basis
            
            m = size(self.V,2);
            md = size(U,2);
            [~,I] = find(SIZE == m);
            TDouble = zeros(SIZE); % if md=1, array of dimension mxmxmd are computed as matrices mxm, but we want a tensor mxmx1 as a final result
            T = tenzeros(SIZE);
            Elements = self.Mesh.Elements;
            V = self.V;                %#ok<*PROPLC>
            
            % parsing element weights
            [elementWeights,inputs] = self.parse_inputs(varargin{:});
            
            % extracting elements with nonzero weights
            elementSet = find(elementWeights);
            
            % Computing element level contributions

            for j = elementSet
                thisElement = Elements(j).Object;
                index = thisElement.iDOFs;          
                Ve = V(index,:); 
                Ue = U(index,:);
                if strcmpi(mode,'ELP') % toggle Element-Level projection
                    Te = thisElement.(elementMethodName)(inputs{1}(j,:),inputs{2}, inputs{3}, inputs{4});
                    Ter = einsum('iI,ijk,jJ,kK->IJK',Ve,Te,Ve,Ue);
                    TDouble = TDouble + Ter;
                else
                      msg = 'Only the ELP mode is currently supported for hydrodynamic forces';
                      error(msg)
                end
            end

            if md==1
                T(:,:,1) = TDouble;    
            else
                T = tensor(TDouble);
            end 
           
            
            % Next possible improvement: store the tensor as a sparse
            % tensor
            %[subs, T] = sparsify(T,[],sumDIMS);
            %T = sptensor(subs, T, SIZE);
        end

        function [T] = tensor4_skin_PROM(self,elementMethodName,U,SIZE,sumDIMS,mode,varargin)
            % This function assembles a generic finite element vector from
            % its element level counterpart.
            % elementMethodName is a string input containing the name of
            % the method that returns the element level vector Fe.
            % For this to work, a method named elementMethodName which
            % returns the appropriate vector must be defined for all
            % element types in the FE Mesh.            
            % NOTE: in this function, we reduce all the dimensions with 
            % the same reduction basis
            
            m = size(self.V,2);
            md = size(U,2);
            [~,I] = find(SIZE == m);
            TDouble = zeros(SIZE); % if md=1, array of dimension mxmxmd are computed as matrices mxm, but we want a tensor mxmx1 as a final result
            T = tenzeros(SIZE);
            Elements = self.Mesh.Elements;
            V = self.V;                %#ok<*PROPLC>
            
            % parsing element weights
            [elementWeights,inputs] = self.parse_inputs(varargin{:});
            
            % extracting elements with nonzero weights
            elementSet = find(elementWeights);
            
            % Computing element level contributions

            for j = elementSet
                thisElement = Elements(j).Object;
                index = thisElement.iDOFs;          
                Ve = V(index,:); 
                Ue = U(index,:);
                if strcmpi(mode,'ELP') % toggle Element-Level projection
                    Te = thisElement.(elementMethodName)(inputs{1}(j,:),inputs{2}, inputs{3}, inputs{4});
                    Ter = einsum('iI,ijkl,jJ,kK,lL->IJKL',Ve,Te,Ve,Ve,Ue);
                    TDouble = TDouble + Ter;
                else
                      msg = 'Only the ELP mode is currently supported for hydrodynamic forces';
                      error(msg)
                end
            end

            if md==1
                T(:,:,:,1) = TDouble;    
            else
                T = tensor(TDouble);
            end 
           
            
            % Next possible improvement: store the tensor as a sparse
            % tensor
            %[subs, T] = sparsify(T,[],sumDIMS);
            %T = sptensor(subs, T, SIZE);
        end

    end
end
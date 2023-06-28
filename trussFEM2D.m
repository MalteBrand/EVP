classdef trussFEM2D

    methods (Static)

        function [u,f,Ke,K] = solve(k,b,EAs,BCs,loads)
            % k: knots, matrix: [x, y; ...
            % b: bars, matrix: [knot 1, knot 2; ...
            % BCs: boundary conditions, matrix: [knot index, DOF; ...
            % EAs: bar stiffness per bar, vector
            % loads: loads, matrix: [knot index, DOF, magnitude; ...
                
            nDOF = length(k)*2; % total degrees of freedom before BCs
            nob = length(b); % number of bars
            
            % element stiffness matrices and assembly of global
            Ke = zeros(4,4,nob);
%             K = zeros(nDOF,nDOF);
            indi=zeros(4,4,nob);
            indj=zeros(4,4,nob);
            for i=1:nob
                Ke(:,:,i) = trussFEM2D.barElement2D(k(b(i,1),1),k(b(i,1),2),k(b(i,2),1),k(b(i,2),2),EAs(i));
                I1 = b(i,1)*2-1; 
                I2 = b(i,2)*2-1;
                indii = [I1:I1+1, I2:I2+1];
                indi(:,:,i)=repmat(indii,4,1);
                indj(:,:,i)=repmat(indii,4,1)';
%                 K(indii,indii) = K(indii,indii) + Ke(:,:,i);
            end
            K=sparse(indi(:),indj(:),Ke(:),nDOF,nDOF);
            % load vector
            f = zeros(nDOF,1);
            for i=1:size(loads,1)
                f( loads(i,1)*2-2+loads(i,2) ) =f( loads(i,1)*2-2+loads(i,2) ) + loads(i,3);
            end
            % boundary conditions
            inBC=zeros(length(BCs(:,1)),1);
            for i=1:length(BCs(:,1))
                inBC(i) = BCs(i,1)*2-2+BCs(i,2);
            end
            inBC = sort(inBC);
            inFree = 1:nDOF;
            inFree(inBC) = [];
            K(inBC,:)=[];
            K(:,inBC)=[];
            ffree=f;
            ffree(inBC) = [];
            % solve
            uFree=K\ffree;
            % add BCs to u
            u=zeros(nDOF,1);
            u(inFree)=uFree;
        end

        function [k,b] = truss_preProcess2D(n,m,dx,dy)
            % structure - n x m
            % nodes / knots
            k = zeros(n*m,2);
            ii=1;
            for i=0:1:n-1
                for j=0:1:m-1
                    k(ii,:) = [i*dx,j*dy];
                    ii=ii+1;
                end
            end
            % bars
            b = zeros(n*m,2);
            ii=1;
            for j=1:n
                for i=1:m
                    if i<m
                        b(ii,:) = [i+(j-1)*m,(i+1)+(j-1)*m];
                        ii=ii+1;
                    end
                    if j<n
                        b(ii,:) = [i+(j-1)*m,i+j*m];
                        ii=ii+1;
                    end
                    if i<m && j<n
                        b(ii,:) = [i+(j-1)*m,i+1+j*m];
                        ii=ii+1;
                        b(ii,:) = [i+1+(j-1)*m,i+j*m];
                        ii=ii+1;
                    end
                end
            end
        end

        function plotTruss2D(k,b,t,fn,u,scaleU)
            % k: node coordinates
            % b: element connection
            % t: thickness of a beam
            % fn: figure index
            % u: deformation
            % scaleU: scale factor for deformation
            scaleT=1e3;
            if nargin>4
                if nargin<6
                    scaleU = ( max(max(k))-min(min(k)) ) / max(abs(u)) *0.1;
                end
                k = k+scaleU*[u(1:2:end), u(2:2:end)];
            end
            myFig=figure(fn); myFig.Color='white';cla; hold on;
            for i=1:length(b)
                plot([k(b(i,1),1),k(b(i,2),1)],[k(b(i,1),2),k(b(i,2),2)], ...
                    'Color','black','LineWidth',t(i)*scaleT,'Marker','o',  'MarkerFaceColor','white');
            end
            axis equal;
            axis off;
        end

        function [Keg,L] = barElement2D(x1,y1,x2,y2,EA)
            % x1, y1: coordinates of node 1
            % x2, y2: coordinates of node 2
            % EA: axial stiffness
            alp = atan((y2-y1)./(x2-x1));
            L = sqrt((x2-x1).^2 + (y2-y1).^2);
            T = [cos(alp) sin(alp) 0 0; 0 0 cos(alp) sin(alp)];
            Kl = EA/L*[1 -1; -1 1];
            Keg = T'*Kl*T;
        end

        function [m,dm]=mass(k,b,As)
            % k: node coordinates
            % b: element connection
            % A: cross section Area
            x1=k(b(:,1),1);
            y1=k(b(:,1),2);
            x2=k(b(:,2),1);
            y2=k(b(:,2),2);
            L = sqrt((x2-x1).^2 + (y2-y1).^2);
            m=sum(As.*L);
            dm=L;
        end

        function [c,dc,u]=compliance(k,b,EAs,BCs,loads)
            % k: node coordinates
            % b: element connection
            % EA: bar stiffness per bar
            % BCs: boundary conditions, matrix: [knot index, DOF; ...
            % loads: loads, matrix: [knot index, DOF, magnitude; ...
            [u,f,Ke]=trussFEM2D.solve(k,b,EAs,BCs,loads);
            c=u'*f;
            nob=size(b,1);
            dc=zeros(nob,1);
            for i=1:nob
                dKe=Ke(:,:,i)/EAs(i);
                I1 = b(i,1)*2-1; 
                I2 = b(i,2)*2-1;
                ind = [I1:I1+1, I2:I2+1];
                ue=u(ind);
                dc(i) = -ue'*dKe*ue;
            end
        end

    end

end

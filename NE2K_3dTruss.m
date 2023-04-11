function K=NE2K_3dTruss(N, E, e, A)
% given Node matrix N, element matrix E, elasicity modulus e, and
% cross-section area A, generate the stiffness matrix K
Nnode=size(N,1);% number of nodes
K=zeros(Nnode*3);% initial global stiffness matrix
Ne=size(E,1);% number of elements
for ie=1:Ne
    xi=N(E(ie,1),1);
    yi=N(E(ie,1),2);
    zi=N(E(ie,1),3);
    xj=N(E(ie,2),1);
    yj=N(E(ie,2),2);
    zj=N(E(ie,2),3);
    L=sqrt((xj-xi)^2+(yj-yi)^2+(zj-zi)^2); % length of the element
    Cx=(xj-xi)/L;
    Cy=(yj-yi)/L;
    Cz=(zj-zi)/L;
    w = [Cx*Cx Cx*Cy Cx*Cz ; Cy*Cx Cy*Cy Cy*Cz ; Cz*Cx Cz*Cy Cz*Cz]; 
    k = e*A(ie)/L*[w -w ; -w w]; % the element stiffness matrix, to be assembled in to K
    K=SpaceTrussAssemble(K, k, E(ie,1), E(ie,2));
end
end

function y = SpaceTrussAssemble(K,k,i,j) 
%SpaceTrussAssemble   This function assembles the element stiffness 
%                     matrix k of the space truss element with nodes 
%                     i and j into the global stiffness matrix K. 
%                     This function returns the global stiffness   
%                     matrix K after the element stiffness matrix   
%                     k is assembled. 
K(3*i-2,3*i-2) = K(3*i-2,3*i-2) + k(1,1); 
K(3*i-2,3*i-1) = K(3*i-2,3*i-1) + k(1,2); 
K(3*i-2,3*i) = K(3*i-2,3*i) + k(1,3); 
K(3*i-2,3*j-2) = K(3*i-2,3*j-2) + k(1,4); 
K(3*i-2,3*j-1) = K(3*i-2,3*j-1) + k(1,5); 
K(3*i-2,3*j) = K(3*i-2,3*j) + k(1,6); 
K(3*i-1,3*i-2) = K(3*i-1,3*i-2) + k(2,1); 
K(3*i-1,3*i-1) = K(3*i-1,3*i-1) + k(2,2); 
K(3*i-1,3*i) = K(3*i-1,3*i) + k(2,3); 
K(3*i-1,3*j-2) = K(3*i-1,3*j-2) + k(2,4); 
K(3*i-1,3*j-1) = K(3*i-1,3*j-1) + k(2,5); 
K(3*i-1,3*j) = K(3*i-1,3*j) + k(2,6); 
K(3*i,3*i-2) = K(3*i,3*i-2) + k(3,1); 
K(3*i,3*i-1) = K(3*i,3*i-1) + k(3,2); 
K(3*i,3*i) = K(3*i,3*i) + k(3,3); 
K(3*i,3*j-2) = K(3*i,3*j-2) + k(3,4); 
K(3*i,3*j-1) = K(3*i,3*j-1) + k(3,5); 
K(3*i,3*j) = K(3*i,3*j) + k(3,6); 
K(3*j-2,3*i-2) = K(3*j-2,3*i-2) + k(4,1); 
K(3*j-2,3*i-1) = K(3*j-2,3*i-1) + k(4,2); 
K(3*j-2,3*i) = K(3*j-2,3*i) + k(4,3); 
K(3*j-2,3*j-2) = K(3*j-2,3*j-2) + k(4,4); 
K(3*j-2,3*j-1) = K(3*j-2,3*j-1) + k(4,5); 
K(3*j-2,3*j) = K(3*j-2,3*j) + k(4,6); 
K(3*j-1,3*i-2) = K(3*j-1,3*i-2) + k(5,1); 
K(3*j-1,3*i-1) = K(3*j-1,3*i-1) + k(5,2); 
K(3*j-1,3*i) = K(3*j-1,3*i) + k(5,3); 
K(3*j-1,3*j-2) = K(3*j-1,3*j-2) + k(5,4); 
K(3*j-1,3*j-1) = K(3*j-1,3*j-1) + k(5,5); 
K(3*j-1,3*j) = K(3*j-1,3*j) + k(5,6); 
K(3*j,3*i-2) = K(3*j,3*i-2) + k(6,1); 
K(3*j,3*i-1) = K(3*j,3*i-1) + k(6,2); 
K(3*j,3*i) = K(3*j,3*i) + k(6,3); 
K(3*j,3*j-2) = K(3*j,3*j-2) + k(6,4); 
K(3*j,3*j-1) = K(3*j,3*j-1) + k(6,5); 
K(3*j,3*j) = K(3*j,3*j) + k(6,6); 
y = K; 
end
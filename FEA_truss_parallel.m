function g=FEA_truss_parallel(e,P)
A1=2; 
A2=1.2;
A3=0.6; 
L=120;  % L=R/2
A=A3*ones(1, 52);% cross-sectional area
A([1:8, 29:36])=A1;
A(9:16)=A2;
S1=sqrt(3/2); 
S2=sqrt(2); 
S3=sqrt(3);

Node=...% the coordinates of all nodes
    [0,0,2*L; 
     -L,0,S3*L;
     0,-L,S3*L;
     L,0,S3*L;  
     0,L,S3*L;
     -S3*L,0,L;
     -S1*L,-S1*L,L;
     0,-S3*L,L; 
     S1*L,-S1*L,L;
     S3*L,0,L;  
     S1*L,S1*L,L;
     0,S3*L,L; 
     -S1*L,S1*L,L;
     -2*L,0,0; 
     -S2*L,-S2*L,0;
     0,-2*L,0;  
     S2*L,-S2*L,0;
     2*L,0,0;   
     S2*L,S2*L,0;
     0,2*L,0;   
     -S2*L,S2*L,0];

Element=...% nodes of all elements
    [1,2;   %1
    1,3;   %2
    1,4;   %3
    1,5; %4
     2,6;   %5
     3,8;   %6
     4,10;  %7
     5,12;%8
     2,7;   %9
     3,7;   %10
     3,9;   %11
     4,9;%12
     4,11;  %13
     5,11;  %14
     5,13;  %15
     2,13; %16
     2,3;   %17
     3,4;   %18
     4,5;   %19
     5,2;  %20
     6,7;   %21
     7,8;   %22
     8,9;   %23
     9,10;%24
     10,11; %25
     11,12; %26
     12,13; %27
     13,6;%28
     6,14;  %29
     7,15;  %30
     8,16;  %31
     9,17;%32
     10,18; %33
     11,19; %34
     12,20; %35
     13,21;%36
     6,15;  %37
     8, 15;  %38
     8,17;  %39
     10,17;%40
     10,19; %41
     12,19; %42
     12,21; %43
     6, 21;%44
     7, 14;  %45
     7,16;  %46
     9,16;  %47
     9,18;%48
     11,18; %49
     11,20; %50
     13,20; %51
     13,14]; %52

 n=length(e);% how many FEM calculations
 g=zeros(size(e));
 P=-P;% all the load follows -Z direction , so minus the magnitude.
 
 for j=1:n
K=NE2K_3dTruss(Node, Element, e(j), A); % assemble the global stiffness matrix
% add displacement constraints
K_index=true(size(K));
K_index(1:39,1:39)=false;
K(K_index)=0;

for i=40:63
K(i,i)=1;
end

% define the load vector/matrix

f=zeros(63,1);% load vector

f(3*(1:13))=P(j,:);
u=K\f;
g(j)=u(3);% the z displacement of Node 1 
 end






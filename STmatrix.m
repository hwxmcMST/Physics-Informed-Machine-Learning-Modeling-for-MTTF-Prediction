function M=STmatrix(z1,z2,z3)
% this function fetch one element of each column of z to form a row, and
% then return a matrix M. If z is a 2 by 3 matrix, then M is a 2^3 by 3
% matrix.
if nargin==3
[n1,m1]=size(z1);
[n2,m2]=size(z2);
[n3,m3]=size(z3);
Nrow=n1*n2*n3;
q=1;
M=zeros(Nrow,m1+m2+m3);
for i=1:n1
    for j=1:n2
        for k=1:n3
            M(q,:)=[z1(i,:),z2(j,:),z3(k,:)];
            q=q+1;
        end
    end
end
end
if nargin==2
    [n1,m1]=size(z1);
    [n2,m2]=size(z2);
    Nrow=n1*n2;
    q=1;
    M=zeros(Nrow,m1+m2);
    for j=1:n1
        for k=1:n2
            M(q,:)=[z1(j,:),z2(k,:)];
            q=q+1;
        end
    end
end
    
            
            
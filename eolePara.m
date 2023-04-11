function [eigv,eigf,cmat]=eolePara(Nterm,node,corr_fun)
%Generate eivenvalue and eigen funciton and correlation mat for random
%fields
S=size(node);
N=S(1);% discretized into N nodes
cmat=ones(N);


if S(2)==1 % one-dimensional random field
for i=1:N
    for j=1:N
        cmat(i,j)=corr_fun([node(i),node(j)]);
    end
end


else % at leat two-dimensional random field
for i=1:N
    for j=1:N
        if j>i
          for k=1:S(2)
          cmat(i,j)=cmat(i,j)*corr_fun{k}([node(i,k),node(j,k)]);                 %([node(i),node(j)]);
          end
        end
        if j<i
            cmat(i,j)=cmat(j,i);
        end
    end
end
    
end

A=eye(N);
[V,lumbda]=eigs(cmat,A,Nterm);
[eigv,idx]  = sort(diag(lumbda),'descend');   
eigf=V(:,idx);
end
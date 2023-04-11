function R=BC2R(B,C,ns)

[V,lumbda]=eig(C);
[eigv,idx]  = sort(diag(lumbda),'descend');
eigv(eigv<1e-4*max(eigv))=[];% remove those small eigen values.
Nterm=length(eigv);
eigf=V(:,idx(1:Nterm));
U=normrnd(0,1,ns,Nterm);
X=repmat(B',ns,1);
X=X+U*sqrt(diag(eigv))*eigf';
n=length(B);
R=B;
Xmin=X(:,1);
for i=1:n
    if i==1
        Xmin=X(:,1);
    else
        Xmin=min([Xmin,X(:,i)],[],2);
    end
    R(i)=sum(Xmin>0)/ns;
end


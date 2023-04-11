function y=eoleOutput(U,z,znode,eigv,eigf,cmat,mu,sig)
% if z is specified, output the realization of ranodm process only at specific z,
% otherwise output the entire realization of random process at specific
% realization U of the random inputs.

[NN,~]=size(znode);

if ~isempty(z)
    [N,~]=size(z);
    y=ones(N,1)*mu;
    Nterm=length(eigv);
    for i=1:N
        T=znode-repmat(z(i,:),NN,1);
        T=T.^2;
        TT=sum(T,2);
        [~,k]=min(TT);% locate the index of z
        for j=1:Nterm
            y(i)=y(i)+sig*U(i,j)/sqrt(eigv(j))*((eigf(:,j))'*(cmat(:,k)));
        end
    end
else
    Ns=length(U(:,1));
%     Nterm=length(eigv);
    Nt=length(cmat(:,1));
    y=mu*ones(Ns,Nt);
    
%     for k=1:Nt
%         for j=1:Nterm
%             y(:,k)=y(:,k)+sig*U(:,j)/sqrt(eigv(j))*((eigf(:,j))'*(cmat(:,k)));
%         end
%     end
    y=y+sig*U*sqrt(diag(eigv))*eigf';
end










% The following codes are used in the first version. To add a
% multi-dimensional random field, as required by the reviewer, I modified
% the codes to the above version. Xinpeng Wei, 8/26/2019
% if ~isempty(t)
%     N=length(t);
%     y=ones(N,1)*mu;
%     Nterm=length(eigv);
%     for i=1:N
%         [~,k,~]=intersect(tnode,t(i));
%         if isempty(k)
%         distance=(tnode-t(i)).^2;
%         [~,k]=min(distance);
%         end
%         for j=1:Nterm
%             y(i)=y(i)+sig*U(i,j)/sqrt(eigv(j))*((eigf(:,j))'*(cmat(:,k)));
%         end
%     end
% else
%     Ns=length(U(:,1));
%     Nterm=length(eigv);
%     Nt=length(cmat(:,1));
%     y=mu*ones(Ns,Nt);
%     for k=1:Nt
%         for j=1:Nterm
%             y(:,k)=y(:,k)+sig*U(:,j)/sqrt(eigv(j))*((eigf(:,j))'*(cmat(:,k)));
%         end
%     end
% end
%     

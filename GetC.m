function C=GetC(alpha,Ts, F_corr)
% calculate the correlation matrix of the equivalent Gaussian process
n=length(Ts);
alpha1=alpha(:,1:4);
alpha2=alpha(:,end);% the last column corresponds to random process F(t)
c=alpha1*alpha1';
C=zeros(size(c));

for i=1:n
    for j=1:n
          rho=F_corr([Ts(i),Ts(j)]);
          C(i,j)=c(i,j)+alpha2(i)*rho*alpha2(j);
    end
end
        
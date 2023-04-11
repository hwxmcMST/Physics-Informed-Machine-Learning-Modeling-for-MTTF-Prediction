function [Rt,Ncall]=FormMCS(p)
% Convert limit-state function into Gaussian process and then sampling the
% Gaussian process to obtain time dependent reliability Rt
%This function is specifically design for Example 2. Do not expect it to
% work for other examples!
NX=5;% in Example 2, 4 random variables + 1 random process
alpha=zeros(p.Nt,NX);
beta=zeros(p.Nt,1);
global functioncall;
functioncall=0;
for i=1:p.Nt
    [alpha(i,:),beta(i)]=FORM(p, p.Ts(i)); % FORM at every time instant
end
C=GetC(alpha,p.Ts,p.corr_F);
Rt=BC2R(beta,C,p.Nmcs);
Ncall=functioncall;
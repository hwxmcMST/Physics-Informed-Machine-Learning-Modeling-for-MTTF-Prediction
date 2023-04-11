function [W,X]=getWX_SG(k,disttype, distpara)
[NX,~]=size(distpara);
[p,W]=sg('GQN',NX,k);
% p=p/sqrt(2);
X=u2x(p,disttype, distpara);



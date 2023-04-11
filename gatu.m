function [g,gradient_u]=gatu(model,u,t,disttype,distpara)
%Calculate g and its 1st & 2nd derivatives analytically at point u (in u space)
%Xiaoping Du, 12/22/2006
n=length(u);
x=u2x(u,disttype,distpara);
g=feval(model,t,x);
step(1:n)=min(1.0e-3,max(abs(u/100),1.0e-5));

    for i=1:n
         temp=u(i); u(i)=u(i)+step(i); 
         x=u2x(u,disttype,distpara);
         gtemp=feval(model,t,x);
         gradient_u(i)=(gtemp-g)/step(i);
         u(i)=temp; 
    end 


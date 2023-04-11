function y= exp3_1(x,t,F)

r1=x(:,1);
r2=x(:,2);
r3=x(:,3);
r4=x(:,4);
rho1 = x(:,5);
rho2 = x(:,6);
rho3 = x(:,7);
rho4 = x(:,8);



rt=(r1.*exp(-0.05*t)).^2;

y=pi.*rt.*rho1-0.707*F;

end
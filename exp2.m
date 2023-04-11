function y=exp2(x,t,F)
density=7.85e4;
k=5e-4;
L=5;
x1=x(:,1);
x2=x(:,2);
x3=x(:,3);
q=x(:,4);
% at=x2-2*k*t; 
% bt=x3-2*k*t;
at=x2.*exp(-0.02*t);
bt=x3.*exp(-0.02*t);
% min(at)
% min(bt)
% if min(at)<0||min(bt)<0
%     error('negative size!');
% end
y=(at.*bt.^2.*x1)/4-(q*L^2/8+F*L/4+(density*x2.*x3*L^2)/8);

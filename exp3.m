function r= exp3(x, t)
E=x(:,1);% elasticity modulus
P=x(:,2:end);% all loads
P=[P,zeros(size(P,1),8)];% Nodes 6~13 are subjected to nothing, in this example
P(:,1)=P(:,1).*exp(0.02*t);% only the mean value of F1 varies with time.
r=1.3+FEA_truss_parallel(E,P);% the displacement is along -z direction, so the returned value of FEA_truss_parallel is negative.
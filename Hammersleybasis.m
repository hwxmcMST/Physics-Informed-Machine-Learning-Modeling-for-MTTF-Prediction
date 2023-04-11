function y=Hammersleybasis(d,k)
%d is the dimension of variable; k is the index of sample
Hbasis=primes(10000);
Hb=Hbasis(d-1);
Hbnew=Hb;knew=k;ynew=0;
for i=1:1:100000
if knew>0
rk=mod(knew,Hb);
ynew=ynew+rk/Hbnew;
knew=floor(knew/Hb);
Hbnew=Hbnew*Hb;
else
break    
end
end
y=ynew;

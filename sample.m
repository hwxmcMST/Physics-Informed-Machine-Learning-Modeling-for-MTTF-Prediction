function x=sample(numsimu,disttype,distpara,option)
[nrandom,column]=size(distpara);

if option==1
    u= lhsgen(numsimu,nrandom);  %Agus' code
end
if option==2
    u=lhsdesign(numsimu,nrandom,'criterion','correlation'); %Matlab code
end
if option==3
    u=rand(numsimu,nrandom); %MCS
end

if option==4 %optimized LHS, onlu for Examples 1 and 2
    if numsimu==100
        fid = fopen('L100-80.txt','r');
    end
    if numsimu==250
        fid = fopen('L250-80.txt','r');
    end
    if numsimu==1000
        fid = fopen('L1000.txt','r');
    end
    for i=1:numsimu
        for j=1:nrandom    % !!!!!!!!!!!!!! 
            s = fscanf(fid,'%s',1);
            u(i,j)=(str2num(s) - 0.5)/numsimu;
            if abs(u(i,j)-1.0)<1.e-6
                u(i,j)=1.0-1.e-6;
            end
            if u(i,j)<=1.0e-6
                u(i,j)=1.e-6;
            end
        end
    end
fclose(fid);
end
% fclose(fid);

for i=1:nrandom
    z=u(:,i);
    p(1:column)=distpara(i,1:column);
    type=disttype(i,:);
    if strcmp(deblank(type),'Normal') | strcmp(deblank(type),'norm')
        x(:,i)=norminv(z,p(1),p(2));
    end
    if strcmp(deblank(type),'Lognormal') | strcmp(deblank(type),'logn')
        b=(log((p(2)/p(1))^2+1))^0.5;
        a=log(p(1))-0.5*b^2;
        x(:,i)=logninv(z,a,b);
    end
    if strcmp(deblank(type),'Weibull') | strcmp(deblank(type),'weib')
        x(:,i)=wblinv(z,p(1),p(2));
    end
    if strcmp(deblank(type),'Exponential') | strcmp(deblank(type),'exp')
        x(:,i)= expinv(z,p(1));
    end
    if strcmp(deblank(type),'Uniform') | strcmp(deblank(type),'unif')
        x(:,i)=unifinv(z,p(1),p(2));
    end
    if strcmp(deblank(type),'Extreme1') | strcmp(deblank(type),'ext1')
        x(:,i)=evinv(z,p(1),p(2));
     end
 end
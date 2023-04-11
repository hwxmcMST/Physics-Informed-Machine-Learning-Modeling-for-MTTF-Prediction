function para=exp3_in()


ExampleID=3;
ContainRP=1;% whether this example contains random process or random field
model=['exp3_1';'exp3_2';'exp3_3';'exp3_4'];%;'exp3_3'

% input random variables
distpara=[2.5,0.05;
          2.5,0.05;
          2.3,0.05;
          2.3,0.05;
          36E3,3E3;
          36E3,3E3;
          36E3,3E3;
          36E3,3E3;
  ];% parameters of input variables
disttype=char('norm','norm','norm','norm','norm','norm','norm','norm');

% discretization of time
tbds=[1,10];% boundary of t
Nt=100;% t is discretized into 50 points
Ts=linspace(tbds(1),tbds(2),Nt)';% discretization points of t

% input random process
corl_F=6;% correlation length of F(t)
corr_F=@(x)exp(-(x(1)-x(2))^2/corl_F^2);% correlation function
muF=500e3;% mean value of F(t)
sigF=50e3;% standard deviation of F(t)
NtermF=6;% use NtermF random variables to represent the process T(t)
[eigvF,eigfF,cmatF]=eolePara(NtermF,Ts,corr_F);% EOLE

% add the random variables from EOLE
distpara=[distpara;repmat([0,1],NtermF,1)];
disttype=char(disttype,repmat(char('norm'),NtermF,1));
NX=length(disttype(:,1));

Nmcs=1e6;% sample size for MCS
Nmcs0=1e3;
Uth=2;
Stdthd=5e-3;
AddEachLoop=1e3;% increase the Nmcs0 by at most 1e2 each loop
Xs=sample(Nmcs0,disttype,distpara,3);% the MCS population of X, i.e. x^MCS in the paper.
Nxi=20;% initial sample points of X
u=(Hammersley(NX+1,Nxi))';% Hammersley to get samples in hypercube space
Ti=u(:,1)*(tbds(2)-tbds(1))+tbds(1);% initial samples of time t
Xi=u2x(norminv(u(:,2:end)), disttype, distpara);% map the initial sample from hypercube to physical space
Xi1=Xi(:,1:(end-NtermF));% original input variables part
Xi2=Xi(:,(end-NtermF+1):end);% input random process part
Fi=eoleOutput(Xi2,Ti,Ts,eigvF,eigfF,cmatF,muF,sigF);

rho_th=0.95;
for i=1:4
Yi(:,i)=feval(model(i,:), Xi1, Ti,Fi);
end




para=struct('u',u,'model',model,'disttype',disttype,'distpara',distpara,'Xs',Xs,'Ts',Ts,'Xi',Xi,'Ti',Ti,'tbds',tbds,...
    'Nmcs',Nmcs,'ContainRP',ContainRP,'Nxi',Nxi,'Yi',Yi,'Uth',Uth, 'ExampleID', ExampleID,...
    'eigvF',eigvF,'eigfF',eigfF,'cmatF',cmatF,'NtermF',NtermF,'Nt',Nt,'muF',muF,'sigF',sigF,'Fi',Fi,'rho_th',rho_th,...
    'Nmcs0', Nmcs0, 'Stdthd',Stdthd,'AddEachLoop',AddEachLoop,'corr_F',corr_F);
end

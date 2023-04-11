function para=exp1_in()
ExampleID=1;
ContainRP=0;% whether this example contains random process or random field
model='exp1';
distpara=[0,1];% parameters of input variables
disttype='unif';% distributio types of input variables
NX=size(distpara,1);% number of input random variable
% note that this number is adaptively increased until convergence.
Nmcs=1e5;% sample size for MCS
Nmcs0=1e3;
Uth=2;
Stdthd=5e-3;
AddEachLoop=1e3;% increase the Nmcs0 by at most 1e2 each loop
Xs=sample(Nmcs0,disttype,distpara,3);
% k=4;% the accuracylevel of sparse grid method
% [W, Xs]=getWX_SG(k,disttype,distpara);
Nxi=5;% initial sample points of X
tbds=[0,40];%the boundaries of t;
Nt=100;%the discretization points of t;
Ts=linspace(tbds(1), tbds(2), Nt);
Ts=Ts';
u=(Hammersley(NX+1,Nxi))';% Hammersley to get samples in hypercube space
Ti=u(:,1)*(tbds(2)-tbds(1))+tbds(1);% initial samples of time t
Xi=u(:,2)*(distpara(2)-distpara(1))+distpara(1);
Yi=feval(model, Xi,Ti);
rho_th=0.95;
para=struct('u',u,'model',model,'disttype',disttype,'distpara',distpara,'Xs',Xs,'Ts',Ts,'Xi',Xi,'Ti',Ti,'tbds',tbds,...
    'Nmcs',Nmcs,'ContainRP',ContainRP,'AddEachLoop',AddEachLoop,'Nxi',Nxi,'Yi',Yi,'Uth',Uth, 'ExampleID', ExampleID,...
    'Nmcs0', Nmcs0,'Stdthd',Stdthd,'rho_th',rho_th,'Nt',Nt);


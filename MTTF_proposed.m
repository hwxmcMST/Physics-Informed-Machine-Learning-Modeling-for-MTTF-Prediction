function [opt,para]=MTTF_proposed(para)

if para.ContainRP==0 % 0 means this example does not contain input random field
    NX=length(para.Xi(1,:));% dimension of X
    NR=NX+1;% number of input variables including t
    regpoly_fun='regpoly0';
    Corr_fun='corrgauss';
    theta = 1*ones(1,NR); lob = 1e-3*ones(1,NR); upb = 1e3*ones(1,NR);
    Nx=length(para.Xs(:,1));% sample size of X
    fttf=zeros(Nx,1);% first time to failure
    sig_fttf=zeros(Nx,1);
    X_train=para.Xi;
    T_train=para.Ti;
    In_train=[X_train,T_train];
    Y_train=para.Yi;
    AddIter=0;
    MTTF=[];
    AddSample=0;
    SampleSize=para.Nmcs0;
    Gamma=[];
    x_den=x2den(para.Xs,para.disttype,para.distpara);
    smooth=1;
while 1
    while 1
        [Y_hat,~]=dacefit([X_train,T_train],Y_train,regpoly_fun, Corr_fun, theta, lob, upb);
        sites=[repmat(para.distpara(:,1),para.Nt,1),para.Ts];
        
        ERP_samples=ERPD(Y_hat,sites,para.Nterm_ERP,para.Ns_ERP);
        
        for k=1:Nx
            Xs_k=para.Xs(k,:);
            Inputk=[repmat(Xs_k,para.Nt,1),para.Ts];
            [muk,sigk]=predictors(Inputk,Y_hat);
            [~,ia,~]=intersect(Inputk,In_train,'rows');
            sigk(ia)=0;
            sigk=sqrt(sigk);
            ERPk=ERP_samples.*repmat(sigk',para.Ns_ERP,1)+repmat(muk',para.Ns_ERP,1);% shift the ERP_samples by specific mean and varance.
            [fttf(k),sig_fttf(k),~]=ERP2FTTF(para.Ts,ERPk,smooth);        
        end
        CurrentMTTF=mean(fttf)
        MTTF=[MTTF;CurrentMTTF]
        [~,X_new_index]=max(x_den.*sig_fttf);%density*standard deviation  is the learning fucntion to determine x
        Input_new=[repmat(para.Xs(X_new_index, :),para.Nt,1),para.Ts];
        [mu_x_new, sig_x_new]=predictors(Input_new,Y_hat);
         [~,ia,~]=intersect(Input_new,In_train,'rows');
        sig_x_new(ia)=0;       
        sig_x_new=sqrt(sig_x_new);
        ERP_x_new=ERP_samples.*repmat(sig_x_new',para.Ns_ERP,1)+repmat(mu_x_new',para.Ns_ERP,1);% shift the ERP_samples by specific mean and varance.
        [~,~,fttf_samples]=ERP2FTTF(para.Ts,ERP_x_new,smooth); 
        [t_new, Pmin]=get_t_new(para.Ts,fttf_samples);
       Pmin
       S=mean(sig_fttf)/CurrentMTTF
        if S<1e-2% stopping criterion
            Pmin
            break
        end
        Xnext=para.Xs(X_new_index,:);
        Ynext=feval(para.model,Xnext,t_new);
        X_train=[X_train;Xnext];
        Y_train=[Y_train;Ynext];
        T_train=[T_train;t_new];
        In_train=[In_train;Xnext,t_new];
        AddIter=AddIter+1
    end
   fttf_std=std(fttf);
   gamma=fttf_std/CurrentMTTF/sqrt(Nx);
   Gamma=[Gamma;gamma];
   if gamma<para.Stdthd
       break;
   end
   AddSample=AddSample+1
   Nx_add=(fttf_std/CurrentMTTF/para.Stdthd)^2-Nx;% need to add Nx_add samples
   Nx_add=round(Nx_add)
   if Nx_add>para.AddEachLoop
   Nx_add=para.AddEachLoop;
   end
   if Nx_add<1
       break;
   end
   Xs_add=sample(Nx_add,para.disttype,para.distpara,3);% increase the population of X
   para.Xs=[para.Xs;Xs_add];% update the population of X
%    para.Nmcs=para.Nmcs+Nx_add;
   SampleSize=[SampleSize;SampleSize(end)+Nx_add];
   Nx=Nx+Nx_add;
   fttf=zeros(Nx,1);
   sig_fttf=zeros(Nx,1);
   x_den=x2den(para.Xs,para.disttype,para.distpara);
end
  opt=struct('fttf',fttf,'CandidatePoints',Nx,'InitialPoints',...
      para.Nxi,'AddedPoints',AddIter,'TrainingPoints',para.Nxi+AddIter,'X_train',...
      X_train,'Y_train',Y_train,'T_train',T_train,'MTTF',MTTF);
end


if para.ExampleID==2 
    NR=6;% five random variables/process and 1 time variable
    NX1=4; % four input random variables
    regpoly_fun='regpoly0';
    Corr_fun='corrgauss';
    theta = 1*ones(1,NR); lob = 1e-3*ones(1,NR); upb = 1e3*ones(1,NR);
    Nx=length(para.Xs(:,1));% sample size of X
    fttf=zeros(Nx,1);% first time to failure
    sig_fttf=zeros(Nx,1); 
    X_train=para.Xi;% including the random process part
    X1_train=para.Xi(:,1:NX1);% does not include the random process part
    T_train=para.Ti;
    F_train=para.Fi;
    Y_train=para.Yi;
    In_train=[X1_train,T_train,F_train];
    AddIter=0;
    MTTF=[];
    AddSample=0;
    SampleSize=para.Nmcs0;
    Gamma=[]; 
    x_den=x2den(para.Xs,para.disttype,para.distpara);
    smooth=0;
while 1        
    while 1
        [Y_hat,~]=dacefit([X1_train,T_train, F_train],Y_train,regpoly_fun, Corr_fun, theta, lob, upb);
%         Theta_t=Y_hat.theta(end-1);
%          ERP_samples=EpistemicRandomProcessDiscretization(Theta_t, para.tbds(2), para.Nt, para.Nterm_ERP, para.Ns_ERP);
        Fk_mean=eoleOutput(zeros(1,para.NtermF),[],[],para.eigvF,para.eigfF,para.cmatF,para.muF,para.sigF);
        sites=[repmat(para.distpara(1:NX1,1)',para.Nt,1),para.Ts,Fk_mean'];
        ERP_samples=ERPD(Y_hat,sites,para.Nterm_ERP,para.Ns_ERP);
        for k=1:Nx
            Xs_k1=para.Xs(k,1:NX1);
            Xs_k2=para.Xs(k,NX1+1:end);
            Fk=eoleOutput(Xs_k2,[],[],para.eigvF,para.eigfF,para.cmatF,para.muF,para.sigF);            
            Inputk=[repmat(Xs_k1,para.Nt,1),para.Ts,Fk'];
            [muk,sigk]=predictors(Inputk,Y_hat);
            [~,ia,~]=intersect(Inputk,In_train,'rows');
            sigk(ia)=0;
            sigk=sqrt(sigk);
            ERPk=ERP_samples.*repmat(sigk',para.Ns_ERP,1)+repmat(muk',para.Ns_ERP,1);% shift the ERP_samples by specific mean and varance.
            [fttf(k),sig_fttf(k),~]=ERP2FTTF(para.Ts,ERPk,smooth); 
        end
         CurrentMTTF=mean(fttf)
        MTTF=[MTTF;CurrentMTTF];
        [~,X_new_index]=max(x_den.*sig_fttf);%density*standard deviation  is the learning fucntion to determine x
        Xs_new1=para.Xs(X_new_index,1:NX1);
        Xs_new2=para.Xs(X_new_index,NX1+1:end);
        Fk_new=eoleOutput(Xs_new2,[],[],para.eigvF,para.eigfF,para.cmatF,para.muF,para.sigF);
        Input_new=[repmat(Xs_new1,para.Nt,1),para.Ts,Fk_new'];
        [mu_x_new, sig_x_new]=predictors(Input_new,Y_hat);
        [~,ia,~]=intersect(Input_new,In_train,'rows');
        sig_x_new(ia)=0;
        sig_x_new=sqrt(sig_x_new);
        ERP_x_new=ERP_samples.*repmat(sig_x_new',para.Ns_ERP,1)+repmat(mu_x_new',para.Ns_ERP,1);% shift the ERP_samples by specific mean and varance.
        [~,~,fttf_samples]=ERP2FTTF(para.Ts,ERP_x_new,smooth); 
        [t_new, Pmin]=get_t_new(para.Ts,fttf_samples);
        Pmin
      S=mean(sig_fttf)/CurrentMTTF
%       S2=norm(sig_fttf)/length(sig_fttf)/CurrentMTTF
         if S<1e-2% stopping criterion
            break
        end
        
        Xnext=para.Xs(X_new_index,:);
        X1next=Xnext(1:NX1);
        X2next=Xnext(NX1+1:end);
        Fnext=eoleOutput(X2next,t_new,para.Ts,para.eigvF,para.eigfF,para.cmatF,para.muF,para.sigF);
        Ynext=feval(para.model,X1next,t_new,Fnext);
        X_train=[X_train;Xnext];
        X1_train=[X1_train;X1next];
        Y_train=[Y_train;Ynext];
        T_train=[T_train;t_new];
        F_train=[F_train;Fnext];
        In_train=[In_train;X1next,t_new,Fnext];
        AddIter=AddIter+1
    end
    
   fttf_std=std(fttf);
   gamma=fttf_std/CurrentMTTF/sqrt(Nx);
   Gamma=[Gamma;gamma];
   if gamma<para.Stdthd
       break;
   end
   AddSample=AddSample+1
   Nx_add=(fttf_std/CurrentMTTF/para.Stdthd)^2-Nx;% need to add Nx_add samples
   Nx_add=round(Nx_add)
   if Nx_add>para.AddEachLoop
   Nx_add=para.AddEachLoop;
   end
   if Nx_add<1
       break;
   end
   Xs_add=sample(Nx_add,para.disttype,para.distpara,3);% increase the population of X
   para.Xs=[para.Xs;Xs_add];% update the population of X
   SampleSize=[SampleSize;SampleSize(end)+Nx_add];
   Nx=Nx+Nx_add;
   fttf=zeros(Nx,1);
  x_den=x2den(para.Xs,para.disttype,para.distpara);
end
  opt=struct('fttf',fttf,'CandidatePoints',Nx,'InitialPoints',...
      para.Nxi,'AddedPoints',AddIter,'TrainingPoints',para.Nxi+AddIter,'X_train',...
      X_train,'Y_train',Y_train,'T_train',T_train,'MTTF',MTTF);
end

function [opt,para]=KrigMttf_rho(para)



if para.ExampleID==3
    NR=4;% five random variables/process and 1 time variable
    NX1=8; % four input random variables
    regpoly_fun='regpoly0';
    Corr_fun='corrgauss';
    theta = {1*ones(1,NR),1*ones(1,NR),1*ones(1,NR),1*ones(1,NR)}; 
    lob = {1e-3*ones(1,NR),1e-3*ones(1,NR),1e-3*ones(1,NR),1e-3*ones(1,NR)};
    upb = {1e3*ones(1,NR),1e3*ones(1,NR),1e3*ones(1,NR),1e3*ones(1,NR)};
    Nx=length(para.Xs(:,1));% sample size of X
    fttf=zeros(Nx,4);% first time to failure
    LocalMinU=zeros(Nx,2);% Min U learning funcion wrt t, with x frozen
    c1x=para.Xi(:,[1,5]); c2x=para.Xi(:,[2,6]); c3x=para.Xi(:,[3,7]);c4x=para.Xi(:,[4,8]);
    X_train={c1x,c2x,c3x,c4x};% including the random process part
    t1=para.Ti;t2=para.Ti;t3=para.Ti;t4=para.Ti;
    T_train={t1,t2,t3,t4};
    f1=para.Fi; f2=para.Fi;f3=para.Fi;f4=para.Fi;
    F_train={f1,f2,f3,f4};
    y1=para.Yi(:,1); y2=para.Yi(:,2); y3=para.Yi(:,3);y4=para.Yi(:,4);
    Y_train={y1,y2,y3,y4};
    AddIter=0;
    MTTF=[];
    AddSample=0;
    SampleSize=para.Nmcs0;
    Gamma=[]; 
while 1        
    while 1
       for j=1:4
              [Y_hat,~]=dacefit([X_train{j},T_train{j},F_train{j}],Y_train{j},regpoly_fun, Corr_fun, theta{j}, lob{j}, upb{j});
            for k=1:Nx
                Fk=eoleOutput(para.Xs(k,NX1+1:end),[],[],para.eigvF,para.eigfF,para.cmatF,para.muF,para.sigF);
                if j==1
                    [Utemp(:,j),fttf(k,j)]=GetU(para.Xs(k,[1,5]),para.Ts,Fk',Y_hat);
                    [LocalMinU(k,2,j),LocalMinU(k,1,j)]=min(Utemp(:,j)); 
                elseif j==2
                     [Utemp(:,j),fttf(k,j)]=GetU(para.Xs(k,[2,6]),para.Ts,Fk',Y_hat);
                    [LocalMinU(k,2,j),LocalMinU(k,1,j)]=min(Utemp(:,j)); 
                elseif j==3 
                     [Utemp(:,j),fttf(k,j)]=GetU(para.Xs(k,[3,7]),para.Ts,Fk',Y_hat);
                    [LocalMinU(k,2,j),LocalMinU(k,1,j)]=min(Utemp(:,j)); 
                else
                    [Utemp(:,j),fttf(k,j)]=GetU(para.Xs(k,[4,8]),para.Ts,Fk',Y_hat);
                    [LocalMinU(k,2,j),LocalMinU(k,1,j)]=min(Utemp(:,j)); 
                end
            end        
       end
        fttfsys=min([max(fttf(:,1:2),[],2),max(fttf(:,3:4),[],2)],[],2);
        LocalMinU1=LocalMinU(:,:,1);LocalMinU2=LocalMinU(:,:,2);LocalMinU3=LocalMinU(:,:,3);LocalMinU4=LocalMinU(:,:,4);
        localU={LocalMinU1,LocalMinU2,LocalMinU3,LocalMinU4};
        LocalMinU_temp1 = [LocalMinU1(:,2),LocalMinU2(:,2)];
        LocalMinU_temp2 = [LocalMinU3(:,2),LocalMinU4(:,2)];
        [LocalMinUU1(:,2),LocalMinUU1(:,1)]=max(LocalMinU_temp1,[],2);
        [LocalMinUU2(:,2),LocalMinUU2(:,1)]=max(LocalMinU_temp2,[],2);
        ComU = {LocalMinUU1,LocalMinUU2};

        [LocalMinUU(:,2),LocalMinUU(:,1)]=min([LocalMinUU1(:,2),LocalMinUU2(:,2)],[],2);
%         fttf
        CurrentMTTF=mean(fttfsys);
        MTTF=[MTTF;CurrentMTTF];
        [GlobalMinU,Xindex]=min(LocalMinUU(:,2));
        GlobalMinU
        if GlobalMinU>para.Uth 
            break
        end
      
        Sys_index=LocalMinUU(Xindex,1);
%         ComU{Sys_index}
        Ctemp_index=ComU{Sys_index}(Xindex,1);
        switch  Sys_index
            case 1
                if  Ctemp_index==1
                    C_index=1;
                else
                    C_index=2;
                end
            otherwise
                if Ctemp_index==1
                    C_index=3;
                else
                    C_index=4;
                end
        end
        Xnext=para.Xs(Xindex,:);
        Tindex = localU{C_index}(Xindex,1);
        Tnext=para.Ts(Tindex);
%         X1next=Xnext(1:NX1);
        X2next=Xnext(NX1+1:end);
        Fnext=eoleOutput(X2next,Tnext,para.Ts,para.eigvF,para.eigfF,para.cmatF,para.muF,para.sigF);
        Ynext=feval(para.model(C_index,:),Xnext,Tnext,Fnext);
       
        if C_index==1
            col=[1,5];
        elseif C_index==2
            col=[2,6];
        elseif C_index==3
            col=[3,7];
        else
            col=[4,8];
        end
         X_train{C_index}=[X_train{C_index};Xnext(col)];
%         X1_train{C_index} = [X1_train{C_index};X1next];
        Y_train{C_index}=[Y_train{C_index};Ynext];
        T_train{C_index}=[T_train{C_index};Tnext];
        F_train{C_index} = [F_train{C_index};Fnext];
        AddIter=AddIter+1
    end
    
   fttf_std=std(fttfsys);
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
   clear fttf LocalMinUU LocalMinUU1 LocalMinUU2;
%    fttf=zeros(Nx,1);
%    LocalMinU=zeros(Nx,2);
end
  opt=struct('fttf',fttf,'CandidatePoints',Nx,'InitialPoints',...
      para.Nxi,'AddedPoints',AddIter,'TrainingPoints',para.Nxi+AddIter,'X_train',...
      X_train,'Y_train',Y_train,'T_train',T_train,'MTTF',MTTF);
end



function FTTF=MCS(para)

if para.ExampleID==3
Ts=para.Ts;
Xs=sample(para.Nmcs,para.disttype,para.distpara,3);
Xs1=Xs(:,1:(end-para.NtermF));
Xs2=Xs(:,(end-para.NtermF+1):end);
FTTF=zeros(para.Nmcs,1);
Fs=eoleOutput(Xs2,[],[],para.eigvF,para.eigfF,para.cmatF,para.muF,para.sigF);
for k=1:para.Nmcs
     for j=1:4    
         YY(:,j)=feval( para.model(j,:),repmat(Xs1(k,:),para.Nt,1), Ts, Fs(k,:)');    
     end
     fttf_index1=find(YY(:,1)<=0,1);
     fttf_index2=find(YY(:,2)<=0,1);
     fttf_index3=find(YY(:,3)<=0,1);
     fttf_index4=find(YY(:,4)<=0,1);

     if isempty(fttf_index1)
         fttf1=Ts(end);
     else
         fttf1=Ts(fttf_index1);
     end

     if isempty(fttf_index2)
         fttf2 = Ts(end);
     else
         fttf2=Ts(fttf_index2);
     end

     if isempty(fttf_index3)
         fttf3 = Ts(end);
     else
         fttf3=Ts(fttf_index3);
     end

     if isempty(fttf_index4)
         fttf4 = Ts(end);
     else
         fttf4=Ts(fttf_index4);
     end

     FTTF(k)=min(max(fttf1,fttf2),max(fttf3,fttf4));
end
end    

if para.ExampleID==4
Ts=para.Ts;
Xs=sample(para.Nmcs,para.disttype,para.distpara,3);
Xs1=Xs(:,1);
Xs2=Xs(:,2:end);
FTTF=zeros(para.Nmcs,1);
end 
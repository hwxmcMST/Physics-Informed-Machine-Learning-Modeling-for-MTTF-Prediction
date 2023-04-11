function [alpha, beta]=FORM(p,t)

D=5;
u0=zeros(1,D);
u=u0;
beta=norm(u);
iter_max=100;
conv_flag=0;
distpara=[p.distpara(1:4,:);p.muF,p.sigF];
disttype=char(p.disttype(1:4,:),'norm');
i=0;

while conv_flag==0
    [g,grad_g]=gatu(p.model_form,u,t,disttype,distpara); 
    beta0=beta;
    beta=beta+g/norm(grad_g);
    alpha=grad_g/norm(grad_g);
    u=-beta*alpha;
    % Check convergence
      if abs(beta-beta0)/beta0<1e-4
        conv_flag =1;
      end
    if i>=iter_max
        conv_flag =1; 
       disp('Maximum number of iterations is reached for MPP search.');
    end
    i=i+1;
end
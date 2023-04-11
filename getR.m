function R=getR(x_doe, x_try, Ssc, theta, corr_name)
% calculate the correlation matrix of Kriging model
[m,n] = size(x_doe);  % number of design sites and number of dimensions
  [mx,~] = size(x_try);            % number of trial sites and their dimension
  % Normalize trial sites  
  x = (x_try - repmat(Ssc(1,:),mx,1)) ./ repmat(Ssc(2,:),mx,1);

    dx = zeros(mx*m,n);  
    kk = 1:m;
    for  k = 1 : mx
      dx(kk,:) = repmat(x(k,:),m,1) - x_doe;
      kk = kk + m;
    end
%     xs=repmat(x,1,1,m);
%     xs2=permute(xs,[2,3,1]);
%     %clear xs
%     xs3=reshape(xs2,nx,[]);
%     %clear xs2
%     dx=xs3'-repmat(dmodel.S,mx,1);   
     r = feval(corr_name, theta, dx);
    R = reshape(r, m, mx);
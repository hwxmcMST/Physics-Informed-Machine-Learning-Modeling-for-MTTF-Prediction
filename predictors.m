function  [y, mse] = predictors(x, dmodel)%a simple version of predictor
% global Ncall
% Ncall=Ncall+length(x(:,1));
% x=xx;
% dmodel=kriging;
[m,n] = size(dmodel.S);  % number of design sites and number of dimensions
  [mx,nx] = size(x);            % number of trial sites and their dimension
  % Normalize trial sites  
  x = (x - repmat(dmodel.Ssc(1,:),mx,1)) ./ repmat(dmodel.Ssc(2,:),mx,1);
  y = zeros(mx,1);         % initialize result

    dx = zeros(mx*m,n);  kk = 1:m;
%     for  k = 1 : mx
%       dx(kk,:) = repmat(x(k,:),m,1) - dmodel.S;
%       kk = kk + m;
%     end
    xs=repmat(x,1,1,m);
    xs2=permute(xs,[2,3,1]);
    %clear xs
    xs3=reshape(xs2,nx,[]);
    %clear xs2
    dx=xs3'-repmat(dmodel.S,mx,1);   
    
    % Get regression function and correlation
    f = feval(dmodel.regr, x);
    r = feval(dmodel.corr, dmodel.theta, dx);
    r = reshape(r, m, mx);
    
    % Scaled predictor 
    sy = f * dmodel.beta + (dmodel.gamma * r).';
    % Predictor
    y = repmat(dmodel.Ysc(1,:),mx,1) + repmat(dmodel.Ysc(2,:),mx,1) .* sy;
      dC=full(dmodel.C);
      rt = dC\ r;
      u = dmodel.G \ (dmodel.Ft.' * rt - f.');
      mse = repmat(dmodel.sigma2,mx,1) .* repmat((1 + colsum(u.^2) - colsum(rt.^2))',1,1);
   
%>>>>>>>>>>>>>>>>   Auxiliary function  ====================

function  s = colsum(x)
% Columnwise sum of elements in  x
if  size(x,1) == 1,  s = x; 
else
    s = sum(x); 
end
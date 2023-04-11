function [U, fttf]=GetU(x,Ts,Fs,Krig)
% This function does not need a real y_star, instead, we use the estimated
% maximum from Kriging as the y_star.
xx=repmat(x,length(Ts),1);
input=[xx,Ts, Fs];
% size(input);
[mu_Y,cov_Y]=predictors(input,Krig);
U=abs(mu_Y)./sqrt(cov_Y);
         fttf_index=find(mu_Y<=0,1);
         if isempty(fttf_index)% cannot find any negative value of the output, than fttf is set to the largest Ts
               fttf=Ts(end);
         else
               fttf=Ts(fttf_index);
               if fttf_index<length(Ts) 
                   U(fttf_index+1:end)=3;
               end
         end
end
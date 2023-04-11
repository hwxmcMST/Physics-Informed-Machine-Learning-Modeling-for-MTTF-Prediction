function r=RT2MTTF(R,T)
% turn the survival function R and corresponding T into mean time to
% failure
r=0;
for i=1:length(T)-1
    r=r+(R(i+1)+R(i))*(T(i+1)-T(i))/2;
end
clc;
clear;
para=exp2_in();
[Rt,Ncall]=FormMCS(para);
MTTF=RT2MTTF(Rt,para.Ts);
plot(Rt)
clc;
clear;
rng default
para=exp3_in();
Y=MCS(para);
MTTF_MCS=mean(Y)
n=(std(Y)/MTTF_MCS/0.005)^2;
[OPT,para]=KrigMttf_rho(para);% the proposed method
Proposed=OPT(1).MTTF(end)



err=(Proposed-MTTF_MCS)/MTTF_MCS*100;


%% plot results for example 1
if para.ExampleID==1
    t_interval=para.tbds;
    x_interval=para.distpara;
    t=linspace(t_interval(1),t_interval(2),100);
    x=linspace(x_interval(1), x_interval(2),100);
    [X,T]=meshgrid(x,t);
  Y=exp(-0.05*T).*cos(0.25*T+X);
%   figure;
%   mesh(T,X,Y);
%   hold on;
%   MyPlot=plot3(OPT.T_train(1:para.Nxi), OPT.X_train(1:para.Nxi), OPT.Y_train(1:para.Nxi),'ro' );
%   MyPlot2=plot3( OPT.T_train(para.Nxi+1:end), OPT.X_train(para.Nxi+1:end), OPT.Y_train(para.Nxi+1:end),'b+');
%   legend('Actual limit-state function','Initial training points','Added training points');
%   MyPlot.LineWidth=2;
%   MyPlot2.LineWidth=2;
%   hold off;
  figure;
  contour(T,X,Y,'ShowText','on');
  hold on;
  MyPlot3=plot(OPT.T_train(1:para.Nxi), OPT.X_train(1:para.Nxi),'ro');
  MyPlot4=plot(OPT.T_train(para.Nxi+1:end), OPT.X_train(para.Nxi+1:end),'b+');
  MyPlot3.LineWidth=2;
  MyPlot4.LineWidth=2;
  legend(' Actual contours', 'Initial training points', 'Added training points');
  xlabel('{\itt}');
  ylabel('{\itX}');
  set(gca,'FontSize',15,'FontName','TimesNewRoman');
end


disp('%%%%%%%%%%%%%%%%%%%%%%%%%% Result Display%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('MTTF by MCS is :')
MTTF_MCS
disp('MTTF by the proposed method is:')
Proposed
disp('Relative error is (%):')
err
disp('Final sample size of the proposed method is:')
size(para.Xs,1)
disp('Function evaluations of the proposed method is:')
OPT.TrainingPoints
disp('Sample size of MCS is:')
para.Nmcs
disp('Function evaluations of MCS is:')
para.Nmcs*para.Nt

clc;
clear;
x=0:0.01:15;
A=6;
B=4;
y=wblpdf(x,A,B);
myplot=plot(x,y,'k-');
myplot.LineWidth=2;
  xlabel('Product lifetime (years)');
  ylabel('PDF');
  set(gca,'FontSize',15,'FontName','TimesNewRoman');
  
clc;
clear;
t=linspace(0,40,100);
x=linspace(0,1,100);
[X,T]=meshgrid(x,t);
Y=exp(-0.05*T).*cos(0.25*T+X);
figure;
mesh(T,X,Y);
figure;
contour(T,X,Y,'ShowText','on')
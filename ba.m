clc
clear all
close all
f=@(x,y) cos(1.5*(x.^2+y)).*sin(1.5*(x+y.^2));
%f=@(x,y) cos(x.*y);
num_pts=50; 
sigma=0.05;
[x,y]=meshgrid(linspace(-1,1,num_pts),linspace(1,-1,num_pts));
z=f(x,y)+normrnd(0,sigma,size(x));
x=x(:) ; y=y(:); z=z(:);%z=f(x,y)+ normrnd(0,sigma,size(x));
figure(2)
plot3(x,y,z,'r.')
hold on
%number of elements
x_seg=10;%divide x into how many part
y_seg=10;%divide y into how many part
getsurf2(x_seg,y_seg,x,y,z,num_pts);
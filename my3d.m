clc
clear all
close all
%artitrary function
xx=@(u,v) cos((v/2)*pi/2)*cos((u/2)*pi);
yy=@(u,v) cos(((v/3)*pi/2))*sin(u/4)*pi;
zz=@(u,v) sin(v*pi/2);

%noise
sigma=0.1;
num_pts=50;
x_seg=10;%divide x into how many part
y_seg=10;%divide y into how many part
%artitrary points
[u,v]=meshgrid(linspace(-1,1,num_pts),linspace(1,-1,num_pts));
x = xx(u,v)+ normrnd(0,sigma,size(u));
y = yy(u,v)+ normrnd(0,sigma,size(u));
z = zz(u,v)+ normrnd(0,sigma,size(u));
figure(1)
%mesh(x,y,z);
plot3(x,y,z,'r.')
hold on

u=u(:);v=v(:);x=x(:) ; y=y(:); z=z(:);
figure(2)
plot3(x,y,z,'r.','Markersize',3)
hold on
% getsurf2(x_seg,y_seg,u,v,x,num_pts);
% getsurf2(x_seg,y_seg,u,v,y,num_pts);
% getsurf2(x_seg,y_seg,u,v,z,num_pts);
getsurf3d(x_seg,y_seg,u,v,x,y,z,num_pts);
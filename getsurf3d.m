function [] = getsurf3d( x_seg,y_seg,u,v,x,y,z,num_pts)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
u_nodes=linspace(min(u),max(u), x_seg+1)'; % separate from min to max in x elements
v_nodes=linspace(max(v),min(v), y_seg+1)'; % separate from min to max in y elements 
groupx=[u v x];
groupy=[u v y];
groupz=[u v z];
%basis functions 
psi_1=@(x1,x2,y1,y2,x,y)(x-x2).*(y-y2)./((x1-x2).*(y1-y2));
psi_2=@(x1,x2,y1,y2,x,y)-(x-x1).*(y-y2)./((x1-x2).*(y1-y2));
psi_3=@(x1,x2,y1,y2,x,y)(x-x1).*(y-y1)./((x1-x2).*(y1-y2));
psi_4=@(x1,x2,y1,y2,x,y)-(x-x2).*(y-y1)./((x1-x2).*(y1-y2));
%global normal matrix_type
N=zeros((x_seg+1)*(y_seg+1));
nx=zeros((x_seg+1)*(y_seg+1),1);
ny=zeros((x_seg+1)*(y_seg+1),1);
nz=zeros((x_seg+1)*(y_seg+1),1);
ver=zeros(4,1);
verplot=[];
L=[];
for j=1:y_seg
    for i=1:x_seg    
    F1 = find(groupx(:,1)<u_nodes(i+1)&groupx(:,1)>=u_nodes(i));      %find x from Bond,between two x nodes
     B1x = groupx(F1,:);  B1y = groupy(F1,:);B1z = groupz(F1,:);                                             % 
     F2 = find(B1x(:,2)>=v_nodes(j+1)&B1x(:,2)<=v_nodes(j));          %find y from B1,between two y nodes
     B2x = B1x(F2,:);B2y = B1y(F2,:);B2z = B1z(F2,:);
    
    u_ind=B2x(:,1);v_ind=B2x(:,2);x_ind=B2x(:,3);y_ind=B2y(:,3);z_ind=B2z(:,3);
    
     u1=u_nodes(i); %start , match the control point data to the elements
     u2=u_nodes(i+1); %end noded
     v1=v_nodes(j); %start , match the control point data to the elements
     v2=v_nodes(j+1); %end noded
     A=[psi_1(u1,u2,v1,v2,u_ind,v_ind) psi_2(u1,u2,v1,v2,u_ind,v_ind)...
        psi_3(u1,u2,v1,v2,u_ind,v_ind) psi_4(u1,u2,v1,v2,u_ind,v_ind)]; %set up and calculate the design matrix A
     N_element=A'*A;
     nx_element=A'*x_ind;
     ny_element=A'*y_ind;
     nz_element=A'*z_ind;
     %Adding to the global normal equation system
     ver(1)=(j-1)*(x_seg+1)+i;
     ver(2)=ver(1)+1;
     ver(3)=ver(1)+x_seg+2;
     ver(4)=ver(3)-1;
     verplot=[verplot;ver];
     for p=1:4
         nx(ver(p))=nx(ver(p))+nx_element(p);
         ny(ver(p))=ny(ver(p))+ny_element(p);
         nz(ver(p))=nz(ver(p))+nz_element(p);
         for q=1:4
            N(ver(p),ver(q))=N(ver(p),ver(q))+N_element(p,q);           
         end         
     end       
    end    
end

cond(N)
x_nodes=N\nx;
y_nodes=N\ny;
z_nodes=N\nz;
%plot the spline
c=0;

 for j=1:y_seg
    for i=1:x_seg
    c=c+1;
    m=4*(c-1)+1;
%coordinates of the points

U=linspace(u_nodes(i),u_nodes(i+1),num_pts/x_seg);
V=linspace(v_nodes(j),v_nodes(j+1),num_pts/y_seg);
[uu,vv]=meshgrid(U,V);
 u1=u_nodes(i); %start , match the control point data to the elements
     u2=u_nodes(i+1); %end noded
      v1=v_nodes(j); %start , match the control point data to the elements
     v2=v_nodes(j+1); %end noded
xplot=x_nodes(verplot(m))*psi_1(u1,u2,v1,v2,uu,vv)+x_nodes(verplot(m+1))*psi_2(u1,u2,v1,v2,uu,vv)...
    +x_nodes(verplot(m+2))*psi_3(u1,u2,v1,v2,uu,vv)+x_nodes(verplot(m+3))*psi_4(u1,u2,v1,v2,uu,vv);
yplot=y_nodes(verplot(m))*psi_1(u1,u2,v1,v2,uu,vv)+y_nodes(verplot(m+1))*psi_2(u1,u2,v1,v2,uu,vv)...
    +y_nodes(verplot(m+2))*psi_3(u1,u2,v1,v2,uu,vv)+y_nodes(verplot(m+3))*psi_4(u1,u2,v1,v2,uu,vv);
zplot=z_nodes(verplot(m))*psi_1(u1,u2,v1,v2,uu,vv)+z_nodes(verplot(m+1))*psi_2(u1,u2,v1,v2,uu,vv)...
    +z_nodes(verplot(m+2))*psi_3(u1,u2,v1,v2,uu,vv)+z_nodes(verplot(m+3))*psi_4(u1,u2,v1,v2,uu,vv);
figure(2)
mesh(xplot,yplot,zplot)

hold on

    end
 end

end


function [] = getsurf2( x_seg,y_seg,x,y,z,num_pts)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
x_nodes=linspace(min(x),max(x), x_seg+1)'; % separate from min to max in x elements
y_nodes=linspace(max(y),min(y), y_seg+1)'; % separate from min to max in y elements 
group=[x y z];
%basis functions 
psi_1=@(x1,x2,y1,y2,x,y)(x-x2).*(y-y2)./((x1-x2).*(y1-y2));
psi_2=@(x1,x2,y1,y2,x,y)-(x-x1).*(y-y2)./((x1-x2).*(y1-y2));
psi_3=@(x1,x2,y1,y2,x,y)(x-x1).*(y-y1)./((x1-x2).*(y1-y2));
psi_4=@(x1,x2,y1,y2,x,y)-(x-x2).*(y-y1)./((x1-x2).*(y1-y2));
%global normal matrix_type
N=zeros((x_seg+1)*(y_seg+1));
n=zeros((x_seg+1)*(y_seg+1),1);
ver=zeros(4,1);
verplot=[];
L=[];
for j=1:y_seg
    for i=1:x_seg    
    F1 = find(group(:,1)<x_nodes(i+1)&group(:,1)>=x_nodes(i));      %find x from Bond,between two x nodes
     B1 = group(F1,:);                                               % 
     F2 = find(B1(:,2)>=y_nodes(j+1)&B1(:,2)<=y_nodes(j));          %find y from B1,between two y nodes
     B2 = B1(F2,:);
    
    x_ind=B2(:,1);y_ind=B2(:,2);z_ind=B2(:,3);
 
     x1=x_nodes(i); %start , match the control point data to the elements
     x2=x_nodes(i+1); %end noded
      y1=y_nodes(j); %start , match the control point data to the elements
     y2=y_nodes(j+1); %end noded
     A=[psi_1(x1,x2,y1,y2,x_ind,y_ind) psi_2(x1,x2,y1,y2,x_ind,y_ind)...
        psi_3(x1,x2,y1,y2,x_ind,y_ind) psi_4(x1,x2,y1,y2,x_ind,y_ind)]; %set up and calculate the design matrix A
     N_element=A'*A;
     n_element=A'*z_ind;
     L=[L;z_ind];
     %Adding to the global normal equation system
     ver(1)=(j-1)*(x_seg+1)+i;
     ver(2)=ver(1)+1;
     ver(3)=ver(1)+x_seg+2;
     ver(4)=ver(3)-1;
     verplot=[verplot;ver];
     for p=1:4
         n(ver(p))=n(ver(p))+n_element(p);
         for q=1:4
            N(ver(p),ver(q))=N(ver(p),ver(q))+N_element(p,q);           
         end         
     end       
    end    
end

cond(N)
z_nodes=N\n;
%plot the spline
c=0;
L_p=[];
 for j=1:y_seg
    for i=1:x_seg
    c=c+1;
    m=4*(c-1)+1;
%coordinates of the points
X=linspace(x_nodes(i),x_nodes(i+1),num_pts/x_seg);
Y=linspace(y_nodes(j),y_nodes(j+1),num_pts/y_seg);
[x,y]=meshgrid(X,Y);
 x1=x_nodes(i); %start , match the control point data to the elements
     x2=x_nodes(i+1); %end noded
      y1=y_nodes(j); %start , match the control point data to the elements
     y2=y_nodes(j+1); %end noded
zplot=z_nodes(verplot(m))*psi_1(x1,x2,y1,y2,x,y)+z_nodes(verplot(m+1))*psi_2(x1,x2,y1,y2,x,y)...
    +z_nodes(verplot(m+2))*psi_3(x1,x2,y1,y2,x,y)+z_nodes(verplot(m+3))*psi_4(x1,x2,y1,y2,x,y);
figure(2)
mesh(x,y,zplot)
hold on
    end
 end
[px,py]=meshgrid(x_nodes,y_nodes);
px=reshape(px',[(x_seg+1)*(y_seg+1),1]);
py=reshape(py',[(x_seg+1)*(y_seg+1),1]);
plot3(px,py,z_nodes,'k*')
hold on
end


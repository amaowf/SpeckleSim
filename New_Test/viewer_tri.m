
clear all
close all
clc

filename = 'pp_N40_L2.0um_010.txt';
p = dlmread(filename);

filename = 'tt_N40_L2.0um_010.txt';
t = dlmread(filename);
Color=[0.75 0.75 0.75]; 
% Color='g';
% 
% Cut=find (t(4,:)>1);
% t(:,Cut)=[];
t = t';
p=p';
for m= 1:length(t)
    n=t(1:3,m);
    X(1:3,m)=p(1,n)';
    Y(1:3,m)=p(2,n)';
    Z(1:3,m)=p(3,n)';
end
g=fill3(X,Y,Z,Color);
Fsize = 20;    
axis('equal')
xlabel('x (µm)','FontSize', Fsize, 'FontName','Arial');
ylabel('y (µm)','FontSize', Fsize, 'FontName','Arial');
rotate3d on
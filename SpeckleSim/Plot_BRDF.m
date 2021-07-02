
close all
clear all
clc

file_name = 'Result_Efield_tri_BRDF.txt';
a = dlmread(file_name);
m = length(a);
d_theta = (75-(-75))/(m-1);
theta_o = -75:d_theta: 75;
sum_a = a*0;
% for n = 2:2
%     if n<10
%         file_name = strcat(sprintf('Result_Efield_tri_BRDF_00%d.txt', n));
%     else 
%         file_name = strcat(sprintf('Result_Efield_tri_BRDF_0%d.txt', n));
%     end
% b = dlmread(file_name);
% 
% sum_a = sum_a + b;
% 
% plot(theta_o, b)
% hold on
% end 
theta_i = -30;
lambda = 0.6e-6;
sigma_s = 170e-9;
Lc = 5.4e-6;
x = BeckmannKirchhoff(sigma_s, Lc, lambda, theta_o, theta_i);
figure
plot(theta_o, sum_a)
hold on
plot(theta_o, a)
% plot(theta_o, x*0.515e-6)

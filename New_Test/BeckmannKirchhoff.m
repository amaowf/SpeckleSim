function  x = BeckmannKirchhoff(sigma_s, Lc, lambda, theta_o, theta_i)
% clc
% close all
% clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input parameters for Beckmann

% lambda = 0.6e-6; %unit: meter
c0=3.0e+8;
f = c0./lambda*1e+9;
omega = 2*pi*f;
Lc = 1.e-6; %correlation length, micrometer range
% sigma_s = 0.030e-6; % rms of the rough surface

k_0 = 2*pi/lambda;
As = 100e-6^2; %illumination area, mu.2

M =2;
i=0;
% theta_i =-30; %degree
Inte_A = 0;
phi = 0;
%TIS = 0.571;
%sigma_s = lambda*sqrt(0-log(1-TIS))/cosd(theta_i)/(4*pi)
TIS = 1- exp(0-(4*pi*cosd(theta_i)*sigma_s/lambda)^2);
% theta_o = -75:1:75;
mm = length(theta_o);
D(mm) = 0;

% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure
%sum_AA = (dlmread(filenamex));
for theta = theta_o
     i = i+1;

%F = (1+ cosd(theta_i)*cosd(theta)-sind(theta_i)*sind(theta)*cosd(phi))/(cosd(theta_i)*(cosd(theta_i)+cosd(theta)));
F = cosd(theta)/cosd(theta_i);
g = (2*pi*sigma_s/lambda)*(cosd(theta_i)+cosd(theta))^2;
Nue_xy = k_0*sqrt(sind(theta_i)^2-2*sind(theta_i)*sind(theta)*cosd(phi)+sind(theta)^2);

    Inte_A = 0;    
    for m = 1:1:100   
       for n = 1:1:m
          M = M*n; %m!
       end 
    Inte_A = g^m/(m*M).*exp(0-(Nue_xy^2*Lc^2/(4*m)))+Inte_A;   
    M =1;
    end
    D(i) = pi*Lc^2*F^2*exp(-g)/As*Inte_A; % Mean_Scattered_Power, D(rou) of quation (5) 

end

x = D;
% figure ('Unit', 'centimeter', 'Position', [1 1 30 18]); %22 for Transmittance with colorbar, 20 for extinction
% ah =axes('Position',[0.25 0.2 0.42 0.7],'LineWidth',1,'FontSize',22, 'FontName','Arial'); 
% 
% plot(theta_o, D*0.3e+5, 'r', 'linewidth', 2);
% hold on
% plot (theta_sie, AA*0.8e-13, 'b', 'linewidth', 2);
% axis([-80 80 0 0.7])
% set(gca,'XTick',[-80 -40 0  40 80], 'TickDir','out','TickLength',[.02 0], 'FontName','Arial');
% 
% xlabel('\theta_s (°)','FontSize', 22, 'FontName','Arial');
% ylabel('BRDF (a.u.)','FontSize',22, 'FontName','Arial');
% title('E_x at z = 0.5 mm', 'FontSize', 22, 'FontName','Arial');
end %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

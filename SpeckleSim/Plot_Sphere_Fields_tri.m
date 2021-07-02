clear all
close all
clc


Nob_x = 201; 
Nob_y = 201;
D = 4;
material = 'Ag';
x_max = 2.5e-6*D;
y_max = 2.5e-6*D;
dx = 2*x_max/(Nob_x-1);
dy = 2*y_max/(Nob_y-1);

x = -x_max:dx:x_max;
y = -y_max:dy:y_max;

lda = 600; 
k0 = 2*pi/lda;
[xx, zz] = meshgrid(y, x);
E_in = 1*exp(-1j*k0*xx);

Field = 'Efield';
Ix = zeros(Nob_x, Nob_y);
Iy = Ix;
Iz_M = Ix;
Ix_M = Ix;
% Result_Efield_tri_xz_100
aa = dlmread(strcat('Result_', Field, '_tri_xz_088.txt'));
% bb_im = dlmread(strcat('Result_', Field, '_tri_imag.txt'));
for i = 1:Nob_y
    for j = 1:Nob_x        
    Ix(j, i) = abs(aa((i-1)*Nob_x + j, 1) + 1j*aa((i-1)*Nob_x + j, 4)).^2;
    Iz(j, i) = abs(aa((i-1)*Nob_x + j, 3) + 1j*aa((i-1)*Nob_x + j, 6)).^2;
    end 
end

% Field = 'Ex';
% Size = num2str(D);
% aa_M = dlmread(strcat('Mie_Re', Field, '_', material,'_', Size, 'nm_lambda', num2str(lda), 'nm.txt'));
% bb_M = dlmread(strcat('Mie_Im', Field, '_', material,'_', Size, 'nm_lambda', num2str(lda), 'nm.txt'));
% Ix_M = abs(aa_M + 1j*bb_M).^2;
% 
% Field = 'Ez';
% aa_M = dlmread(strcat('Mie_Re', Field, '_', material,'_', Size, 'nm_lambda', num2str(lda), 'nm.txt'));
% bb_M = dlmread(strcat('Mie_Im', Field, '_', material,'_', Size, 'nm_lambda', num2str(lda), 'nm.txt'));
% Iz_M = abs(aa_M + 1j*bb_M).^2;
% 
Fsize = 20;
theta = 0: 2 : 360;
R = D/2*1e-9;
x_c = R*sind(theta);
y_c = R*cosd(theta);

[dummy, z_c] = meshgrid(x_c, y_c);
figure ('Unit', 'centimeter', 'Position', [2 2 20 14]);
ah =axes('Position',[0.15 0.2 0.55 0.65],'LineWidth',3,'FontSize',20, 'FontName','Arial');
imagesc(x*1e+6, y*1e+6, (Ix+Iz))

caxis([0.0 1.3])
colorbar
axis([-7.5 7.5 -7.5 7.5])   
% hold on
% plot(x_c*1e+6, y_c*1e+6, '-w', 'linewidth', 1)

xlabel('z (µm)','FontSize', Fsize, 'FontName','Arial');
ylabel('x (µm)','FontSize', Fsize, 'FontName','Arial');
set(gca,'XTick',[-8, -6, -4 -2 0 2 4, 6, 8],  'TickLength',[.02 0], 'FontSize', Fsize);

title('SIE', 'FontWeight', 'Normal', 'FontSize', Fsize);
% name_out = sprintf(strcat('SIE_EIntensity_', Size, 'nm_',material,'.jpg'));
% print('-djpeg', name_out, '-r300');

% figure ('Unit', 'centimeter', 'Position', [2 2 20 14]);
% ah =axes('Position',[0.15 0.2 0.55 0.65],'LineWidth',3,'FontSize',20, 'FontName','Arial');
% imagesc(x*1e+6, y*1e+6, Ix_M+Iz_M)
% caxis([0.0 10.5])
% 
% colorbar
% axis([-x_max*1e+6 x_max*1e+6 -y_max*1e+6 y_max*1e+6])   

% hold on
% plot(x_c*1e+6, y_c*1e+6, '-w', 'linewidth', 1)
% set(gca,'XTick',[-2 -1 0 1 2],  'TickLength',[.02 0], 'FontSize', Fsize);
% set(gca,'YTick',[-2 -1 0 1 2],  'TickLength',[.02 0], 'FontSize', Fsize);
% 
% xlabel('z (µm)','FontSize', Fsize, 'FontName','Arial');
% ylabel('x (µm)','FontSize', Fsize, 'FontName','Arial');
% title('MIE', 'FontWeight', 'Normal', 'FontSize', Fsize);
% name_out = sprintf(strcat('Mie_EIntensity_', Size, 'nm_',material,'.jpg'));
% print('-djpeg', name_out, '-r300');
% 
% %%%%%%%%%%%%%%%%%%%%%%%
% %Compare the results along a line
% figure ('Unit', 'centimeter', 'Position', [2 2 20 14]);
% ah =axes('Position',[0.15 0.2 0.55 0.65],'LineWidth',3,'FontSize',20, 'FontName','Arial');
% p = floor((Nob_x-1)/2);
% x0 = -x_max + p*dx;
% 
% % Along z-direction
% % I = Ix(:, p) + Iz(:, p)*0;
% % I_M = Ix_M(:, p) + Iz_M(:, p)*0;
% % plot(x*1e+6, I, 'linewidth', 2)
% % hold on
% % plot(x*1e+6, I_M, 'linewidth', 2)
% 
% % %Along x-direction
% I = Ix(p, :) + Iz(p, :);
% I_M = Ix_M(p, :) + Iz_M(p, :);
% plot(x*1e+6, I, 'linewidth', 2)
% hold on
% plot(x*1e+6, I_M, 'linewidth', 2)
% 
% xlabel('z (µm)','FontSize', Fsize, 'FontName','Arial');
% ylabel(strcat('E-field Intensity'),'FontSize', Fsize, 'FontName','Arial');
% set(gca,'XTick',[-1 0 1],  'TickLength',[.02 0], 'FontSize', Fsize);
% % axis([-1, 1, 0.0, 5.4])
% h=legend('MIE', 'SIE');
% legend boxoff;
% set(h, 'Fontsize', 16);
% set(h, 'Location', 'northwest');
% 
% tit = strcat(sprintf('at x = %5.2f ', x0*1e+6), ' µm');
% title(tit, 'FontWeight', 'Normal', 'FontSize', Fsize);
% name_out = sprintf(strcat('Line_EIntensity_', Size, 'nm_Au.jpg'));
% print('-djpeg', name_out, '-r300');
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%
% % %Plot error
% % figure ('Unit', 'centimeter', 'Position', [2 2 20 14]);
% % ah =axes('Position',[0.15 0.2 0.55 0.65],'LineWidth',3,'FontSize',20, 'FontName','Arial');
% % plot(x*1e+6, abs((I-I_M)./I), 'linewidth', 2)
% % axis([-1, 1, 0, 0.5])
% % 
% % xlabel('z (µm)','FontSize', Fsize, 'FontName','Arial');
% % ylabel(strcat('Relative error'),'FontSize', Fsize, 'FontName','Arial');
% % set(gca,'XTick',[-1 0 1],  'TickLength',[.02 0], 'FontSize', Fsize);
% % % h=legend('MIE', 'SIE');
% % % legend boxoff;
% % set(h, 'Fontsize', 16);
% % set(h, 'Location', 'northwest');
% % 
% % tit = strcat(sprintf('at x = %5.2f ', x0*1e+6), ' µm');
% % title(tit, 'FontWeight', 'Normal', 'FontSize', Fsize);
% % name_out = sprintf(strcat('Error_Line_EIntensity_x', Size, 'nm_',material,'.jpg'));
% % print('-djpeg', name_out, '-r300');
% % 
% % 

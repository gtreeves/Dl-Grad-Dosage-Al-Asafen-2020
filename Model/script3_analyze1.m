% script3_analyze1 
% This script analyzes all the data sets that are provided to it 
% (o/p of script2_find) and makes cdf plots and ranges of different 
% parameter sets.

clear
clc
close all

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% CONTROLS
Results_folder  = './';

cdfplots   = 'off';
twocompare = 'on';
otherplots = 'off';

addpath('Functions')
load(strcat(Results_folder,'Found_Params'));
%load(strcat(Results_folder,'Found_Errors'));
load(strcat(Results_folder,'Found_Concs'));
 
Params = Found_Params;
Concs  = Found_Concs;
 %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 

% 
% Extracting all Params
%
KeqD      = 4;
KeqDC     = 1;
m         = 51;           % Number of nuclei in NC14  
stage     = 'interphase';
[~, Ac, ~, Vn, Vc] = nuclearSize(1, 'static', m, stage);
lambdaU   = Params(:,1);
lambdaW   = Params(:,2);
gamma     = Params(:,5);
beta      = Params(:,6);
kappa     = Params(:,8);
gamma_by_beta       = gamma./beta;
Rho  = (lambdaW./(Vn*KeqDC + Vc))./(lambdaU./(Vn*KeqD + Vc));

%
% Calculating max and min ranges
%
range_beta    = [min(beta)    max(beta)];
range_gamma   = [min(gamma)   max(gamma)];
range_lambdaU = [min(lambdaU) max(lambdaU)];
range_lambdaW = [min(lambdaW) max(lambdaW)];
range_kappa   = [min(kappa)   max(kappa)];
range_gamma_by_beta = [min(gamma_by_beta) max(gamma_by_beta)];
range_Rho = [min(Rho) max(Rho)];



%
% Extract - u_x0, u_x1, u_max, u_max_pos
%
u_x0      = Concs(:,1:3);
u_x1      = Concs(:,4:6);
u_max_val = Concs(:,7:9);
u_max_pos = Concs(:,10:12);
len       = size(Concs,1);

%
% Length scale ratio - average
%
r2_1    = u_max_val(:,3)./u_max_val(:,2);
r1_1    = u_max_val(:,1)./u_max_val(:,2);
r2_avg1 = mean(r2_1);
r1_avg1 = mean(r1_1);

% How many of the configurations dl0_4x goes to zero at x = 1
config_zeros = u_x1(:,3) < 1e-10;
temp         = u_x1(config_zeros);

number_dl04X_to_zero   = len - size(temp,1)
out_of = len




% PLOTS
if strcmp(cdfplots,'on')
    
    %figure("paperpositionmode", "auto")
    
    cdfplotDU(lambdaU);
    set(gca,'xscale', 'log','FontSize', 20)
    legend('lambdaU');
    title('lambdaU')
    %pbaspect([1,1,1])
    
    %truesize(gcf,[1000, 1000])
    %print('-r500', '-dpng')
    
    
    set(gca,'xscale', 'log','FontSize', 20)
    cdfplotDU(lambdaW);
    title('lambdaW')
    legend('lambdaW');
    set(gca,'xscale', 'log','FontSize', 20)

    
    set(gca,'xscale', 'log','FontSize', 20)
    cdfplotDU(gamma);
    title('gamma')
    legend('gamma');
    set(gca,'xscale', 'log','FontSize', 20)

    
    set(gca,'xscale', 'log','FontSize', 20)
    cdfplotDU(beta);
    title('beta')
    legend('beta');
    set(gca,'xscale', 'log','FontSize', 20)

    
    set(gca,'xscale', 'log','FontSize', 20)
    cdfplotDU(kappa);
    title('\kappa')
    legend('\kappa');
    xlabel('Value of \kappa')
    ylabel('Fraction of values')
    set(gca,'xscale', 'log','FontSize', 20)
    xticks([10^-2 10^-1 10^0])
    
    set(gca,'xscale', 'log','FontSize', 20)
    cdfplotDU(gamma_by_beta);
    title('gamma/beta')
    legend('gamma by beta');
    set(gca,'xscale', 'log','FontSize', 20)

    
    set(gca,'xscale', 'log','FontSize', 20)
    cdfplotDU(Rho);
    title('Ratio of length scales')
    legend('\rho');
    xlabel('Value of \rho')
    ylabel('Fraction of values')
    set(gca,'xscale', 'log','FontSize', 20)
    xticks([10^-2 10^0 10^2 10^4 10^6])
    xlim([10^-3,10^6])
   
    set(gca,'FontSize', 20)
    cdfplotDU(r2_1);
    title('Amplitude 4X/2X');
    legend('4x/2x');
    set(gca,'xscale', 'log','FontSize', 20)
    
    
    set(gca,'FontSize', 20)
    cdfplotDU(r1_1);
    title('Amplitude 1X/2X');
    legend('1X/2X');
    set(gca,'xscale', 'log','FontSize', 20)
 
end

 

if strcmp(twocompare,'on')
    
    figure("paperpositionmode", "auto")
    %set(gca,'xscale', 'log','yscale','log','FontSize', 20, 'MarkerSize', 10)
    plot(gamma, beta,'.','markersize',20)
    title('beta vs gamma')
    set(gca,'xscale', 'log','yscale','log','FontSize', 20)

    figure("paperpositionmode", "auto")
    set(gca,'xscale', 'log','yscale','log','FontSize', 20)
    plot(lambdaU, lambdaW,'.','markersize',20)
    title('lambdaW vs lambdaU')
    set(gca,'xscale', 'log','yscale','log','FontSize', 20)

    figure("paperpositionmode", "auto")
    set(gca,'xscale', 'log','yscale','log','FontSize', 20)
    plot(kappa, beta,'.','markersize',20)
    title('beta vs kappa')
    set(gca,'xscale', 'log','yscale','log','FontSize', 20)


    figure("paperpositionmode", "auto")
    set(gca,'xscale', 'log','yscale','log','FontSize', 20)
    plot(kappa, gamma,'.','markersize',20)
    title('gamma vs kappa')
    set(gca,'xscale', 'log','yscale','log','FontSize', 20)

    figure("paperpositionmode", "auto")
    set(gca,'xscale', 'log','yscale','log','FontSize', 20)
    plot(kappa, lambdaU,'.','markersize',20)
    title('LambdaU vs kappa')
    set(gca,'xscale', 'log','yscale','log','FontSize', 20)

    figure("paperpositionmode", "auto")
    set(gca,'xscale', 'log','yscale','log','FontSize', 20)
    plot(kappa, lambdaW,'.','markersize',20)
    title('lambdaW vs kappa')
    set(gca,'xscale', 'log','yscale','log','FontSize', 20)
    
    figure("paperpositionmode", "auto")
    set(gca,'xscale', 'log','yscale','log','FontSize', 20)
    plot(kappa, Rho,'.','markersize',20)
    title('length Scale Ratio vs kappa')
    set(gca,'xscale', 'log','yscale','log','FontSize', 20)
    
    figure("paperpositionmode", "auto")
    set(gca,'xscale', 'log','FontSize', 20)
    plot(kappa, r2_1,'.','markersize',20)
    title('Amplitude ratio 4X/2X vs \kappa')
    xlabel("\kappa")
    ylabel("Amplitude ratio 4x/2x")
    set(gca,'xscale', 'log','FontSize', 20)
    xticks([10^-2 10^-1 10^0])
    saveas(gcf,'kappa_4x','tif')
    
    figure("paperpositionmode", "auto")
    set(gca,'xscale', 'log','FontSize', 20)
    plot(kappa, r1_1,'.','markersize',20)
    title('Amplitude ratio 1X/2X vs \kappa')
    xlabel("\kappa")
    ylabel("Amplitude ratio 1x/2x")
    set(gca,'xscale', 'log','FontSize', 20)
    xticks([10^-2 10^-1 10^0])
    saveas(gcf,'kappa_1x','tif')
    
    figure("paperpositionmode", "auto")
    plot(Rho, r1_1,'.','markersize',20)
    title('Amplitude ratio 1X/2X vs \rho')
    xlabel("\rho")
    ylabel("Amplitude ratio 1x/2x")
    set(gca,'xscale', 'log','FontSize', 20)
    xticks([10^-2 10^0 10^2 10^4 10^6])
    xlim([10^-3,10^6])
    saveas(gcf,'rho_1x','tif')
    
    figure("paperpositionmode", "auto")
    set(gca,'xscale', 'log','FontSize', 20)
    plot(Rho, r2_1,'.','markersize',20)
    title('Amplitude ratio 4X/2X vs \rho')
    xlabel("\rho")
    ylabel("Amplitude ratio 4x/2x")
    set(gca,'xscale', 'log','FontSize', 20)
    xticks([10^-2 10^0 10^2 10^4 10^6])
    xlim([10^-3,10^6])
    saveas(gcf,'rho_4x','tif')
end


if strcmp(otherplots,'on')
    
    % Plot dl0_4x/2x and dl0_2x/1x
    xaxis     = [ones(len,1) 2*ones(len,1) 3*ones(len,1)];
    u_max_avg = [mean(u_max_val(:,1)) mean(u_max_val(:,2)) mean(u_max_val(:,3))];
    u_max_std = [std(u_max_val(:,1))  std(u_max_val(:,2))  std(u_max_val(:,3))];
    hold off
    scatter(xaxis(:,1),u_max_val(:,1));
    hold on
    scatter(xaxis(:,2),u_max_val(:,2));
    hold on
    scatter(xaxis(:,3),u_max_val(:,3));
    hold on
    errorbar(xaxis(1,1:3),u_max_avg, u_max_std, 'k', 'LineWidth', 2, 'Capsize', 35);
    hold off
    
    % 3D plot!
    plot3(lambdaW, lambdaU, kappa, 'o')
    xlabel('lambdaW')
    ylabel('lambdaU')
    zlabel('kappa')
    set(gca,'XScale','log', 'YScale','log','ZScale','log')
    hold on
    scatter3(lambdaW, lambdaU, kappa)
    hold off
    grid on
end

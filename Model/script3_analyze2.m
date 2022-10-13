% script3_analyze2
% This script analyzes all the data sets that are provided to it and
% makes cdf plots and ranges of different parameter sets. This script in
% addition provides for making comparisons. Variable v helps make
% these comparisons.

clear
clc
close all

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% CONTROLS
Results_folder    = './';

cdfplots   = 'off';
twocompare = 'on';
otherplots = 'off';

addpath('Functions')
load(strcat(Results_folder,'Found_Params'));
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


% Calculating max and min ranges
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
len       = length(Concs);


% Length scale ratio - average
u_max_val(:,1) = u_x0(:,1);
u_max_val(:,2) = u_x0(:,2);
u_max_val(:,3) = u_x0(:,3);

r2_1    = u_max_val(:,3)./u_max_val(:,2);
r1_1    = u_max_val(:,1)./u_max_val(:,2);
r2_avg1 = mean(r2_1);
r1_avg1 = mean(r1_1);

%
% How many of the configurations dl0_4x goes to zero at x = 1
config_zeros = u_x1(:,3) < 1e-10;
temp         = u_x1(config_zeros);
number_dl04X_to_zero   = len - length(temp)


%
% Testing
%
v =  Rho < 1; 

lambdaU2 = lambdaU(v);
lambdaW2 = lambdaW(v);
gamma2   = gamma(v);
beta2    = beta(v);
kappa2   = kappa(v);
gamma_by_beta2       = gamma2./beta2;
Rho2  = (lambdaW2./(Vn*KeqDC + Vc))./(lambdaU2./(Vn*KeqD + Vc));

lambdaU3 = lambdaU(~v);
lambdaW3 = lambdaW(~v);
gamma3   = gamma(~v);
beta3    = beta(~v);
kappa3   = kappa(~v);
gamma_by_beta3       = gamma3./beta3;
Rho3  = (lambdaW3./(Vn*KeqDC + Vc))./(lambdaU3./(Vn*KeqD + Vc));

number_satisfy_criteria = length(lambdaU2)
out_of = length(lambdaU)


u_max_val2 = u_max_val(v,:);
u_max_val3 = u_max_val(~v,:);

r2_2    = u_max_val2(:,3)./u_max_val2(:,2);
r1_2    = u_max_val2(:,1)./u_max_val2(:,2);
r2_3    = u_max_val3(:,3)./u_max_val3(:,2);
r1_3    = u_max_val3(:,1)./u_max_val3(:,2);
r2_avg2 = mean(r2_2);
r1_avg2 = mean(r1_2);
r2_avg3 = mean(r2_3);
r1_avg3 = mean(r1_3);

%
% Plots
%
if strcmp(cdfplots,'on')
    
    cdfplotDU3(lambdaU, lambdaU2, lambdaU3);
    title('lambdaU')
    set(gca,'xscale', 'log')
    legend('All', 'Satisfied', 'Not Satisfied');
   
    cdfplotDU3(lambdaW, lambdaW2, lambdaW3);
    title('lambdaW')
    set(gca,'xscale', 'log')
    legend('All', 'Satisfied', 'Not Satisfied');

    cdfplotDU3(gamma, gamma2, gamma3);
    title('gamma')
    set(gca,'xscale', 'log')
    legend('All', 'Satisfied', 'Not Satisfied');

    cdfplotDU3(beta, beta2, beta3);
    title('beta')
    set(gca,'xscale', 'log')
    legend('All', 'Satisfied', 'Not Satisfied');

    cdfplotDU3(kappa, kappa2, kappa3);
    title('kappa')
    set(gca,'xscale', 'log')
    legend('All', 'Satisfied', 'Not Satisfied');
    
    cdfplotDU3(gamma_by_beta, gamma_by_beta2, gamma_by_beta3);
    title('gamma/beta')
    set(gca,'xscale', 'log')
    legend('All', 'Satisfied', 'Not Satisfied');

    cdfplotDU3(Rho, Rho2, Rho3);
    set(gca,'xscale', 'log')
    title('Ratio of length scales 2')
    legend('All', 'Satisfied', 'Not Satisfied');
    
    cdfplotDU3(r2_1, r2_2, r2_3);
    title('Amplitude 4X/2X');
    legend('All', 'Satisfied', 'Not Satisfied');
    
    cdfplotDU3(r1_1, r1_2, r1_3);
    title('Amplitude 1X/2X');
    legend('All', 'Satisfied', 'Not Satisfied');
    
end

 

if strcmp(twocompare,'on')
    
    figure
    scatter(gamma, beta)
    set(gca,'xscale', 'log','yscale','log')
    title('beta vs gamma')
    hold on
    scatter(gamma2, beta2, 'r')
    hold off

    figure
    plot(lambdaU, lambdaW, 'o')
    set(gca,'xscale', 'log','yscale','log')
    title('lambdaW vs lambdaU')
    hold on
    scatter(lambdaU2, lambdaW2, 'r')
    hold off

    figure
    scatter(kappa, beta)
    set(gca,'xscale', 'log','yscale','log')
    title('beta vs kappa')
    hold on
    scatter(kappa2, beta2, 'r')
    hold off

    figure
    scatter(kappa, gamma)
    set(gca,'xscale', 'log','yscale','log')
    title('gamma vs kappa')
    hold on
    scatter(kappa2, gamma2, 'r')
    hold off
    
    figure
    scatter(kappa, lambdaU)
    set(gca,'xscale', 'log','yscale','log')
    title('LambdaU vs kappa')
    hold on
    scatter(kappa2, lambdaU2, 'r')
    hold off
    
    figure
    scatter(kappa, lambdaW)
    set(gca,'xscale', 'log','yscale','log')
    title('lambdaW vs kappa')
    hold on
    scatter(kappa2, lambdaW2, 'r')
    hold off
    
    figure
    scatter(kappa, Rho)
    set(gca,'xscale', 'log','yscale','log')
    title('length Scale Ratio vs kappa')
    hold on
    scatter(kappa2, Rho2, 'r')
    hold off
    
    figure
    scatter(kappa, r2_1)
    set(gca,'xscale', 'log')
    title('Amplitude 4X/2X vs kappa')
    hold on
    scatter(kappa2, r2_2, 'r')
    hold off
    
    
    figure
    scatter(kappa, r1_1)
    set(gca,'xscale', 'log')
    title('Amplitude 1X/2X vs kappa')
    hold on
    scatter(kappa2, r1_2, 'r')
    hold off
    
    figure
    scatter(Rho, r1_1)
    set(gca,'xscale', 'log')
    title('Amplitude 1X/2X vs Rho')
    hold on
    scatter(Rho2, r1_2, 'r')
    hold off
    
    figure
    scatter(Rho, r2_1)
    set(gca,'xscale', 'log')
    title('Amplitude 4X/2X vs Rho')
    hold on
    scatter(Rho2, r2_2, 'r')
    hold off
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
    scatter3(lambdaW2, lambdaU2, kappa2, 'r')
    hold off
    grid on
end


expplots = true;
if expplots
   figure
   betanew = beta.*Vc*exp(-(0.0001/0.15)^2);
   gammanew = gamma*Vc + gamma*Vn*KeqD*KeqDC;
   bg = betanew./gammanew;
   scatter(Rho,bg)
   xlabel("Rho")
   ylabel("beta_{eff}/gamma_{eff}")
   hold on
   betanew2 = beta2.*Vc*exp(-(0.0001/0.15)^2);
   gammanew2 = gamma2*Vc + gamma2*Vn*KeqD*KeqDC;
   bg2 = betanew2./gammanew2;
   scatter(Rho2, bg2, 'r')
   hold off
   set(gca,'XScale','log')
   ylim([0.4,3])   
end



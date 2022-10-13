% script_run runs parameters to make plots of concentrations
clear
clc
close all

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% CONTROLS
Results_folder = './';
Figs_folder = strcat(Results_folder,'/Figs/');

Testing = 'off';

plotonlydl = "off";
plotboth   = "off"; 
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

addpath('Functions');
if ~exist(Figs_folder, 'dir')
    mkdir(Figs_folder)
end
load(strcat(Results_folder, 'Found_Params'))
load(strcat(Results_folder, 'Found_Concs'))
m       = 51;           % Number of nuclei in NC14
stage   = 'interphase';
[~, Ac, ~, Vn, Vc] = nuclearSize(1, 'static', m, stage);


% 
% Extract Params
%
AllParams = Found_Params;
lambdaU   = AllParams(:,1);
lambdaW   = AllParams(:,2);
KeqD      = 4;
KeqDC     = 1;
gamma     = AllParams(:,5);
beta      = AllParams(:,6);
phi       = 0.15;
kappa     = AllParams(:,8);

gamma_by_beta       = gamma./beta;
Rho  = (lambdaW./(Vn*KeqDC + Vc))./(lambdaU./(Vn*KeqD + Vc));



%
% Initializing variables
%
xspan    = linspace(0,1,m)';
tspan    = linspace(0, 60, 100);
h        = xspan(2) - xspan(1);
reltol   = 1e-7;
abstol   = 1e-7;
fhandle  = @ftns2;
options  = odeset('RelTol',reltol,'AbsTol',abstol,'Jacobian',@jacs2);

%
% Transport Matrix
%
e        = ones(m,1); 
P        = spdiags([e -2*e e],[-1 0 1],m,m);
P(1,2)   = 2; 
P(m,m-1) = 2;


% Testing?
if (strcmp(Testing,'on'))
    v = Rho<1;

    lambdaU2 = lambdaU(v);
    lambdaW2 = lambdaW(v);
    gamma2   = gamma(v);
    beta2    = beta(v);
    kappa2   = kappa(v);
    KeqD2    = repmat(KeqD,length(beta2),1);
    KeqDC2   = repmat(KeqDC,length(beta2),1);
    phi2     = repmat(phi,length(beta2),1);
    Rho2  = (lambdaW2./(Vn*KeqDC + Vc))./(lambdaU2./(Vn*KeqD + Vc));

    AllParams = [lambdaU2, lambdaW2, KeqD2, KeqDC2, gamma2, beta2, phi2, kappa2];
end

n = size(AllParams,1);
error_sna_diff  = zeros(n,2);
error_sogv_diff = zeros(n,2); 
error_sogd_diff = zeros(n,2);
load('Mat/Borders.mat');
fun = @SSE;
options2 = optimset;

for i =2:n      
    
    params = AllParams(i,:);
    [u_space, w_space, u]   = ftn_run(fhandle, tspan, options, params, xspan, P, Vn, Vc, Ac, m);
    [e_sna, e_sogv, e_sogd] = errors_ftn(u_space, xspan);
     error_sna_diff(i,:)  = e_sna;
     error_sogv_diff(i,:) = e_sogv;
     error_sogd_diff(i,:) = e_sogd;
    
    u_space = cell2mat(u_space);
    w_space = cell2mat(w_space);
    [theta_sna]   = fminbnd(fun, 0, 2, options2, xsna,  sigma_xsna,  u_space, xspan);
    [theta_sogv]  = fminbnd(fun, 0, 2, options2, xsogv, sigma_xsogv, u_space, xspan);
    [theta_sogd]  = fminbnd(fun, 0, 2, options2, xsogd, sigma_xsogd, u_space, xspan);

    x1 = [interp1(u_space(:,1), xspan, theta_sna, 'spline');
          interp1(u_space(:,2), xspan, theta_sna, 'spline');
          interp1(u_space(:,3), xspan, theta_sna, 'spline')];


    x2 = [interp1(u_space(:,1), xspan, theta_sogv, 'spline');
          interp1(u_space(:,2), xspan, theta_sogv, 'spline');
          interp1(u_space(:,3), xspan, theta_sogv, 'spline')];

    x3 = [interp1(u_space(:,1), xspan, theta_sogd, 'spline');
          interp1(u_space(:,2), xspan, theta_sogd, 'spline');
          interp1(u_space(:,3), xspan, theta_sogd, 'spline')];
    
      
    % PLOTS   
    if strcmp(plotboth, "on")

        figure("paperpositionmode", "auto")
        plot(xspan, u_space,'--', xspan, w_space, 'LineWidth', 4);        
        xlabel('DV Coordinate', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Concentrations of species', 'FontSize', 12, 'FontWeight', 'bold');
        title('Dl and Dl-Cact distribution', 'FontSize', 14, 'FontWeight', 'bold'); 
        hold on
        plot ([0, x1(1), x1(1)], [theta_sna, theta_sna, 0],'--', 'lineWidth',2 )
        hold on
        plot ([0, x2(1), x2(1)], [theta_sogv, theta_sogv, 0],'--', 'lineWidth',2 )
        hold on
        plot ([0, x3(1), x3(1)], [theta_sogd, theta_sogd, 0],'--', 'lineWidth',2 )
        hold off
        xlim([0,1])
        set(findall(gcf,'-property','FontSize'),'FontSize',20)
        legend('Dl (dl_0 = 0.5)', 'Dl (dl_0 = 1)', 'Dl (dl_0 = 2)','Dl-Cact (dl_0 = 0.5)', 'Dl-Cact (dl_0 = 1)','Dl-Cact (dl_0 = 2)','FontSize',12);
        
        figname = [Figs_folder, 'dl-dl_Cact-', num2str(i),'.jpeg'];
        saveas(gcf,num2str(i),'epsc');
    end
    holder = 1e-3;
    
    if strcmp(plotonlydl, "on")

        figure("paperpositionmode", "auto")
        v = u_space > 10^-3;
        v = double(v);
        u_space = u_space.*v;
        size = 9;

        subplot(1,3,1)
        set(gca,'FontSize',50)
        plot(xspan, u_space, 'LineWidth', 4);
        xlim([0,1])
        legend('Dl (dl_0 = 0.5)', 'Dl (dl_0 = 1)', 'Dl (dl_0 = 2)');
        xlabel('DV Coordinate');
        ylabel('Concentrations of species');
        title('dl distribution');
        set('FontSize', 24)
        
         hold on
         plot ([0, x1(1), x1(1)], [theta_sogv, theta_sogv, holder],'--', 'lineWidth',1.5)
         hold on
         plot ([0, x1(2), x1(2)], [theta_sogv, theta_sogv, holder],'--', 'lineWidth',1.5)
         hold on
         plot ([0, x1(3), x1(3)], [theta_sogv, theta_sogv, holder],'--', 'lineWidth',1.5)
         hold on
         xticks([xsogv])
         xticklabels({})
        
        set(findall(gcf,'-property','FontSize'),'FontSize',20)
        hold off

        %
        subplot(1,3,2)
        semilogy(xspan, u_space, 'LineWidth', 2);
        xlim([0,1])
        legend('Dl (dl0 = 0.5)', 'Dl (dl0 = 1)', 'Dl (dl0 = 2)');
        xlabel('Non-dimensionalized space', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Concentrations of species', 'FontSize', 12, 'FontWeight', 'bold');
        title('Dl distribution', 'FontSize', 14, 'FontWeight', 'bold');
        hold on
        plot ([0, x2(1), x2(1)], [theta_sogv, theta_sogv, holder],'--', 'lineWidth',1 )
        hold on
        plot ([0, x2(2), x2(2)], [theta_sogv, theta_sogv, holder],'--', 'lineWidth',1)
        hold on
        plot ([0, x2(3), x2(3)], [theta_sogv, theta_sogv, holder],'--', 'lineWidth',1)
        xticks([xsogv])
        xticklabels({})
        set(findall(gcf,'-property','FontSize'),'FontSize',size)
        hold off

        subplot(1,3,3)
        semilogy(xspan, u_space, 'LineWidth', 2);
        %ylim([0 0.4])
        xlim([0,1])
        legend('Dl (dl_0 = 0.5)', 'Dl (dl_0 = 1)', 'Dl (dl_0 = 2)');
        xlabel('Non-dimensionalized space', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Concentrations of species', 'FontSize', 12, 'FontWeight', 'bold');
        title('Dl distribution', 'FontSize', 14, 'FontWeight', 'bold');
        %
        hold on
        plot ([0, x3(1), x3(1)], [theta_sogd, theta_sogd, holder],'--', 'lineWidth',0.5 )
        hold on
        plot ([0, x3(2), x3(2)], [theta_sogd, theta_sogd, holder],'--', 'lineWidth',0.5 ) 
        hold on
        plot ([0, x3(3), x3(3)], [theta_sogd, theta_sogd, holder],'--', 'lineWidth',0.5 )
        hold on
        xticks(xsogd)
        xticklabels({})
        %}
        set(findall(gcf,'-property','FontSize'),'FontSize',size)
        hold off
        %}
        figname = [Figs_folder, 'dl-', num2str(i),'.jpeg'];
        %figname = num2str(i);
        saveas(gcf, figname,'epsc');       % Image is the 8th parameter set! 
    end
end

 

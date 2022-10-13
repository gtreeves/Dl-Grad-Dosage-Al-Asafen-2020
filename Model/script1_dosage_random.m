% This script runs the Dorsal dosage model simplified to two equations. and
% stores the parameter values and the corresponding errors on the borders
% of sna and sog. This script picks the parameter values randomly!

clear
clc
close all

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% CONTROLS
%Results_folder = '../Results/Runs2/1/';
Results_folder = './Results/';

N       = 10000;          % Total number of runs 
N_save  = 1000;           % Save results after so many runs

% Ranges of parameters
beta_r    = [10^-3 10^3];
gamma_r   = [10^-3 10^3];
lambdaU_r = [10^-3 10^3];
lambdaW_r = [10^-3 10^3];
kappa_r   = [10^-3 10^3];

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

addpath('Functions');

if ~exist(Results_folder, 'dir')
    mkdir(strcat(Results_folder))
end

% ---------------
% Folder Settings
clk1    = clock;
clk     = [num2str(clk1(1)),  '-', num2str(clk1(2),2),'-', ...
           num2str(clk1(3),2),'_', num2str(clk1(4),2),'-', ...
           num2str(clk1(5),2),'-', num2str(round(clk1(6)),2)];
id      = strcat(clk);
clear clk1 clk;
Mats_folder = strcat(Results_folder,'Mats/', id,'/');
mkdir(Mats_folder); 

%-----------------
% Simulation setup
m       = 51;               % Number of nuclei in NC14
stage   = 'interphase';
[~, Ac, ~, Vn, Vc] = nuclearSize(1, 'static', m, stage);
e        = ones(m,1);       % Contains setup for Transport Matrix
P        = spdiags([e -2*e e],[-1 0 1],m,m);
P(1,2)   = 2; 
P(m,m-1) = 2;
xspan    = linspace(0,1,m)';
tspan    = linspace(0, 60, 500);
h        = xspan(2) - xspan(1);
reltol   = 1e-7;
abstol   = 1e-7;
fhandle  = @ftns2;
options  = odeset('RelTol',reltol,'AbsTol',abstol,'Jacobian',@jacs2);

% ------------------------
% Setting Model parameters
KeqD        = 4;
KeqDC       = 1;
phi         = 0.15;
beta_r      = log10(beta_r);
gamma_r     = log10(gamma_r);
lambdaU_r   = log10(lambdaU_r);
lambdaW_r   = log10(lambdaW_r);
kappa_r     = log10(kappa_r);

%------
% MAIN
AllParams   = zeros(N,8);
AllErrors   = zeros(N,4);
AllConcs    = zeros(N,12);
%negconcs    = zeros(N,1);
%negcount    = 0;
rng shuffle

for i=1:N
    lambdaU = 10^(lambdaU_r(1) + (lambdaU_r(2) - lambdaU_r(1))*rand(1));
    lambdaW = 10^(lambdaW_r(1) + (lambdaW_r(2) - lambdaW_r(1))*rand(1));
    gamma   = 10^(gamma_r(1)   + (gamma_r(2)   - gamma_r(1))*rand(1));
    beta    = 10^(beta_r(1)    + (beta_r(2)    - beta_r(1))*rand(1));
    kappa   = 10^(kappa_r(1)   + (kappa_r(2)   - kappa_r(1))*rand(1));
    params  = [lambdaU, lambdaW, KeqD, KeqDC, gamma, beta, phi, kappa];  
    
    %-----------
    % Run codes
    [u_space, w_space] = ftn_run(fhandle, tspan, options, params, xspan, P, Vn, Vc, Ac, m);
    [error_sna, error_sogv, error_sogd, hard_constraint,u_x0, u_x1, u_max] ... 
        = ftn_errors2(u_space, xspan);

    %-------------
    % Store Values
    AllParams(i,:)   = params;
    AllErrors(i,:)   = [error_sna, error_sogv, error_sogd, hard_constraint];
    AllConcs(i,:)    = [u_x0, u_x1, u_max];
    %{
    if any(cell2mat(u_space(:))<0) || any(cell2mat(w_space(:))<0)
       negconcs(i) = 1;
       negcount    = negcount + 1;
    end
    %}
    %---------------------------------
    % Save results every N_save times.
    if mod(i,N_save) == 0    
        i
        filename = strcat (Mats_folder, 'AllParams', '.mat');
        save(filename, 'AllParams');
        filename = strcat (Mats_folder, 'AllErrors', '.mat');
        save(filename,'AllErrors');
        filename = strcat (Mats_folder, 'AllConcs', '.mat');
        save(filename,'AllConcs');
        if exist("Params_fail","var")
            filename = strcat (Mats_folder, 'Params_fail', '.mat');
            save(filename,'Params_fail');
        end
    end
end














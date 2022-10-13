% script_dosagescaling
%
% This script is set to look at how the dosage sensitivity changes when
% parameters in the empirical (dosage-scaling) model of Dl is used. 
% Those parameters being A,B,phi, and M. 

clear
close all

%
% DV axis
%
nx = 51;
x = linspace(0,1,nx)';
L = 300; % length of DV axis in microns
d = 6; % diameter of one cell

%
% Empirical model parameters
%
phi = 0.15;
m0 = -0.1;
b0 = 0.01 - m0;
alpha00 = 100; % nM
alpha00 = alpha00*0.6; % molc/micron^3

BoverA = 150/400; % from Reeves2012si.pdf, eyeballing Fig S1F,G
bmax = BoverA - 0.7*m0;
nb = 10;
b = linspace(b0,bmax,nb)';

%
% For loop to examine the effect of varying b
%
SxA = zeros(nx,nb);
Shet = zeros(nx,nb);
for i = 1:nb
	b0 = b(i);
	alpha0 = alpha00/(1 + b0); % molc/micron^3
	
	%
	% Dl gradient
	%
	Dl = alpha0*(exp(-x.^2/2/phi^2) + b0 + m0*x);
	
	%
	% Calculate new Dl gradient from change in dosage, then sensitivity
	%
	delt = 1e-2; % differential sensitivity coeff
	Dl2 = Dl*(1 + delt/alpha0);
	x2 = interp1(Dl2,x,Dl); 
	dx = x2 - x;
	SxA(:,i) = alpha0*dx/delt./x;
	
	%
	% Calculate het Dl gradient
	%
	delt = -0.5*alpha0; % het
	Dl2 = Dl*(1 + delt/alpha0);
	x2 = interp1(Dl2,x,Dl);
	dx = x2 - x;
	Shet(:,i) = alpha0*dx/delt./x;

end

%
% Sensitivity coefficient (due to changing dosage)
%
COLR = colorinterpolate(nb,'parula');
figure
for i = 1:nb
	plot(x,SxA(:,i),'Color',COLR{i},'linewidth',2)
	hold on
end
xlim([0 1])
ylim([0 2])
colorbar
set(gcf,'paperpositionmode','auto')
set(gca,'fontsize',24,'xtick',0:0.2:1,'ytick',0:0.5:2)
xlabel('DV coordinate')
ylabel('Sensitivity coeff.')

set(gca,'fontsize',16)








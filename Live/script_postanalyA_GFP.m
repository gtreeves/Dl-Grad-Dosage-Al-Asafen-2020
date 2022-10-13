% script_postanalyA
%
% This script will revise some of the Ahump graphs and try to align them
% all at one point.

clear
close all
colrs = [[0 0.4471 0.7412];[0.8510 0.3255 0.0980];[0.9294 0.6941 0.1255]];
boxplotcolors = [[0 0.4471 0.7412];[0.8510 0.3255 0.0980];[1 1 1];[0.9294 0.6941 0.1255]];
g = [1 2 4];

load 2018-06-14_17-26-32_Soln
load Afit_out
Soln(13) = [];
load A14_avg time A14 Tgast
A_avg = A14;

%
% Empirical, manual estimation of start of nc14, max ampl., and
% squishiness in time.
%
t14 = t14_out;
tgast = tgast_out;
Afit = Afit_out;
delt = (tgast - t14)/Tgast;

%
% Defining the indicies of the different dosages
%
nt = 111;
t = linspace(0,55,nt)';
A = NaN(nt,length(Soln));
G = NaN(length(Soln),1);
M = cell(length(Soln),1);
mednuc = NaN(length(Soln),1);

for i = 1:length(Soln)
	data = Soln(i);
	if strcmp(data.genotype,'dl1x')
		COLR = colrs(3,:);
		G(i) = 1;
	elseif strcmp(data.genotype,'dl2x')
		COLR = colrs(1,:);
		G(i) = 2;
	elseif strcmp(data.genotype,'dl4x')
		COLR = colrs(2,:);
		G(i) = 4;
	else
		COLR = 'k';
	end

end

v1 = G == 1;
v2 = G == 2;
v4 = G == 4;
V = {v1 v2 v4};




%% ========================================================================
% Plot all of the amplitudes as a function of t. Should be aligned to nc14.
% =========================================================================
% {

figure
for i = 1:length(Soln)
	data = Soln(i);
	if strcmp(data.genotype,'dl1x')
		COLR = colrs(3,:);
	elseif strcmp(data.genotype,'dl2x')
		COLR = colrs(1,:);
	elseif strcmp(data.genotype,'dl4x')
		COLR = colrs(2,:);
	else
		COLR = 'k';
	end
	
	plot(data.t-t14(i),data.A,'Color',COLR)
	hold on
	
end
legend('1x','2x','4x')
xlabel('time (min)')
ylabel('Grad. Amplitude (AU)')
set(gca,'Fontsize',16)



%% ========================================================================
% Now, fit the plots to a canonical nc14 curve.
% =========================================================================
% {

figure
for i = 1:length(Soln)
	data = Soln(i);
	if ~isnan(Afit(i))
		plot((data.t-t14(i))/delt(i),data.A/Afit(i))
		hold on
	end
end
plot(time,A_avg,'k','linewidth',2)
xlabel('time (norm)')
ylabel('Grad. Amplitude (norm)')
set(gca,'Fontsize',16)

%}

%% ========================================================================
% An alt to the boxplot representation: point w/errorbars
% =========================================================================
% {

Abar = cellfun(@mean,{Afit(v1) Afit(v2) Afit(v4)});
s = cellfun(@std,{Afit(v1) Afit(v2) Afit(v4)});
% s = s./sqrt(cellfun(@sum,{v1 v2 v4}));

figure
errorbar([1 2 4],Abar,s,'. -','Markersize',16,'color',[0 0 0])
hold on
sig_e = 0.05;
e1 = normrnd(0,sig_e*ones(sum(v1),1));
e2 = normrnd(0,sig_e*ones(sum(v2),1));
e4 = normrnd(0,sig_e*ones(sum(v4),1));
plot(G(v1)+e1,Afit(v1),'.','Color',boxplotcolors(1,:),'Markersize',16)
plot(G(v2)+e2,Afit(v2),'.','Color',boxplotcolors(2,:),'Markersize',16)
plot(G(v4)+e4,Afit(v4),'.','Color',boxplotcolors(4,:),'Markersize',16)


xlim([0.5 4.5])
set(gca,'Fontsize',16)
YLIM = ylim;
ylim([0 1500])
xlabel('dl Copy Number')
ylabel('Max Grad. Amplitude (AU)')


%}

%% ========================================================================
% Given the high degree of variability, run many simulations to get the
% CI's for the live GFP amplitude ratios, using bootstrap
% =========================================================================
% {

N = 5000;
n = cellfun(@sum,V);
I = {randi(n(1),N,n(1)),randi(n(2),N,n(2)),randi(n(3),N,n(3))};

%
% Bootstrap resampling. Using this, we end up with a distribution of
% amplitude ratios (AR).
%
A_boot = zeros(N,length(V));
for i = 1:length(V)
	x = Afit(V{i});
	A_boot(:,i) = mean(x(I{i}),2);
end
	
AR = [A_boot(:,1)./A_boot(:,2) A_boot(:,3)./A_boot(:,2)];

%
% Make histogram
%
figure
X = 0:0.1:3;
histogram(AR(:,1),X,'FaceAlpha',1)
Y12 = histcounts(AR(:,1),X);
hold on
histogram(AR(:,2),X,'FaceAlpha',1)
Y42 = histcounts(AR(:,2),X);

%
% Create normal distributions that match the mean, std, and amplitude of
% the histograms.
%
X1 = X(1:end-1) + diff(X(1:2));
X = linspace(0,3,500);
ARbar = mean(AR);
ARstd = std(AR);
y121 = normpdf(X1,ARbar(1),ARstd(1));
y421 = normpdf(X1,ARbar(2),ARstd(2));
y12 = normpdf(X,ARbar(1),ARstd(1));
y42 = normpdf(X,ARbar(2),ARstd(2));
A12 = sum(Y12.*y121)/sum(y121.^2);
A42 = sum(Y42.*y421)/sum(y421.^2);


%
% Plot the normal distributions on top of the histograms
%
plot(X,A12*y12,'k')%,'color',boxplotcolors(1,:))
plot(X,A42*y42,'k')%'color',boxplotcolors(2,:))

%
% Plot mean and std of distributions on top
%
errorbar(mean(AR(:,1)),100,[],[],std(AR(:,1)),std(AR(:,1)),'k .')
errorbar(mean(AR(:,2)),100,[],[],std(AR(:,2)),std(AR(:,2)),'k .')

legend('1x/2x','4x/2x')
set(gca,'fontsize',24)
xlim([0 3])
xlabel('amplitude ratio')
ylabel('counts')

disp('68% CI:')
disp([mean(AR(:,1))-std(AR(:,1)),mean(AR(:,1))+std(AR(:,1))])
disp([mean(AR(:,2))-std(AR(:,2)),mean(AR(:,2))+std(AR(:,2))])


disp('mean +/- std:')
disp([num2str(mean(AR(:,1))),' +/- ',num2str(std(AR(:,1)))])
disp([num2str(mean(AR(:,2))),' +/- ',num2str(std(AR(:,2)))])



%}


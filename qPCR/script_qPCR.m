% script_qPCR
%
% This script is designed to take all of Allison's qPCR data and make
% weighted means and stdevs, and plot all points w/errorbars and stdev's

clear
close all
gamma = 1;

%
% Load in data
%
[numdata,textdata] = xlsread('Supplementary File 1');

Cq = zeros(9,3,3);
Cq(:,:,1) = numdata(11:19,1:3); % 1x
Cq(:,:,2) = numdata(1:9,1:3); % wt (2x)
Cq(:,:,3) = numdata(21:29,1:3); % 4x

%
% Averaging the technical replicates
%
CQmean1 = zeros(3,3,3);
CQstd1 = zeros(3,3,3);
for o = 1:3
	% each "o" is a different genotype
	for i = 1:3
		% each "i" is a different biological replicate
		CQmean1(i,:,o) = meanDU(Cq(3*(i-1)+(1:3),:,o));
		CQstd1(i,:,o) = stdDU(Cq(3*(i-1)+(1:3),:,o));		
	end
end

%
% Building weighted averages of the biological replicates across different
% days
%
Y = cell(3,1);
S = cell(3,1);
Ybar = zeros(3,1);
Sw = zeros(3,1);
for o = 1:3
	
	%
	% Extract variables
	%
	y = CQmean1(:,:,o); y = y(:);
	s = CQstd1(:,:,o); s = s(:);
	
	%
	% Remove non-numbers, calculate weights
	%
	v = isnan(y) | isinf(y) | isnan(s) | isinf(s);
	y(v) = [];
	s(v) = [];
	Y{o} = y; 
	S{o} = s;
	w = 1./s/sum(1./s);
	
	%
	% Weighted avg
	%
	ybar = sum(w.*y);
	Ybar(o) = ybar;
	
	%
	% Weighted sem
	%
	M = length(y);
	Sw(o) = sqrt(sum(w.*(y - ybar).^2)/(M-1)/mean(w))/sqrt(M);
	
end


%
% Plotting, part 1 (just the means and sem's)
%
figure
x = [1; 2; 4];
c = (2*gamma).^(-Ybar);
actin0 = 1/c(2);
dc = log(2*gamma)*c.*Sw;
errorbar(x,c*actin0,dc*actin0,'k o')
hold on

x1 = linspace(1,4)';
y1 = 0.5*(x1 - 2) + 1;
plot(x1,y1)

xlim([0.5 4.5])
ylim([-0.5 5])
set(gca,'XTick',[1 2 4],'XTicklabel',{'1x','2x','4x'},'Fontsize',24)
print(gcf,'Figs/qpcr_means.eps','-depsc')

set(gca,'Fontsize',16)
xlabel('dl dosage')
ylabel('relative abundance (AU)')
print(gcf,'Figs/qpcr_means.jpg','-djpeg','-r300')


%
% Plotting, part 2 (all data)
%
figure
sig_e = 0.05;
for o = 1:3
	y = Y{o};
	s = S{o};
	
	%
	% Transforming delta-ct into real concentrations.
	%
	c01 = (2*gamma).^(-y);
	dc01 = log(2*gamma)*c01.*s;
	
	e1 = normrnd(0,sig_e*ones(length(y),1));
	errorbar(x(o)+e1,c01*actin0,dc01*actin0,'.')
	hold on
	
end
set(gca,'yscale','log')
hold on

errorbar(x,c*actin0,dc*actin0,'k o')
plot(x1,y1,'k')
xlim([0.5 4.5])
% ylim([0 5])
set(gca,'XTick',[1 2 4],'XTicklabel',{'1x','2x','4x'},'Fontsize',24)
set(gca,'YTick',10.^(-2:2))
print(gcf,'Figs/qpcr_alldata.eps','-depsc')

set(gca,'Fontsize',16)
xlabel('dl dosage')
ylabel('relative abundance (AU)')
print(gcf,'Figs/qpcr_alldata.jpg','-djpeg','-r300')














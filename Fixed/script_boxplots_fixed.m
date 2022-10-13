% script_boxplots_fixed
%
% In this script, we will create boxplots of the Dl gradient widths and
% gene expression boundaries for sna, sog.


clear
close all
gof_gene = 0.8;
gof = 0.7;
boxplotcolors = [[0 0.4471 0.7412];[0.8510 0.3255 0.0980];[0.9294 0.6941 0.1255]];


%% ========================================================================
% Make boxplots of dl grad width in 1x,2x,4x
% =========================================================================
% {

%
% Load sig data (Dl gradient width).
%
load sig_fixed
Sig = {sig1x sigwt sig4x};

%
% We know that dl125 has the wrong shape. So we can do a conversion from 
% the canonical dl125 profile to get an "effective sigma".
%
load refit_dl125 Delt_out gof_out sig1x0
load Dl125_avg
Sig125 = Delt_out*sig1x0;
Sig{1} = Sig125;

%
% Boxplot with the converted 1x
%
boxplotDU(Sig,{'1x','2x','4x'},0.2,boxplotcolors,[],true)
set(gca,'Fontsize',12)
ylabel('Dl gradient width (\sigma)')

%
% Getting slope of log-log plot, which corresponds to the sensitivity
% coefficient of the Dl gradient width with respect to changing dosage.
%
x = [ones(size(Sig125));2*ones(size(Sig{2}));4*ones(size(Sig{3}))];
y = [Sig125;Sig{2};Sig{3}];
[m,b,stats] = linlsq(log2(x),log2(y),[],true);

%}

%% ========================================================================
% Make boxplots of sna,sog gene expression in 1x,2x,4x
% =========================================================================
% {

%
% Load gene expression data
%
load genes_fixed
Wsna = {sD_sna1x sD_snawt sD_sna4x};
SVsog = {sV_sog1x sV_sogwt sV_sog4x};
SDsog = {sD_sog1x sD_sogwt sD_sog4x};

%
% Make boxplots of sog,sna
%
boxplotDU(SVsog,{'1x','2x','4x'},0.2,boxplotcolors,[],true)
set(gca,'Fontsize',12)
ylabel('sog ventral border')

boxplotDU(SDsog,{'1x','2x','4x'},0.2,boxplotcolors,[],true)
set(gca,'Fontsize',12)
ylabel('sog dorsal border')

boxplotDU(Wsna,{'1x','2x','4x'},0.2,boxplotcolors,[],true)
set(gca,'Fontsize',12)
ylabel('sna dorsal border')

%
% Getting slopes of log-log plots
%
x = [ones(size(sV_sog1x));2*ones(size(sV_sogwt));4*ones(size(sV_sog4x))];
y = [sV_sog1x;sV_sogwt;sV_sog4x];
[m,b,stats] = linlsq(log2(x),log2(y),[],true);
x = [ones(size(sD_sog1x));2*ones(size(sD_sogwt));4*ones(size(sD_sog4x))];
y = [sD_sog1x;sD_sogwt;sD_sog4x];
[m,b,stats] = linlsq(log2(x),log2(y),[],true);
x = [ones(size(sD_sna1x));2*ones(size(sD_snawt));4*ones(size(sD_sna4x))];
y = [sD_sna1x;sD_snawt;sD_sna4x];
[m,b,stats] = linlsq(log2(x),log2(y),[],true);


%}

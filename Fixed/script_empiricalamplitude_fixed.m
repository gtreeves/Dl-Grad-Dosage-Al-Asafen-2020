% script_boxplots_fixed
%
% In this script, we will find the amplitudes that, when applied to the
% "canonical" (roughly: average) Dl gradient profiles from 1x, 2x, and 4x
% embryos, result in the most robust possible gene expression boundaries.
%
clear
close all


%% ========================================================================
% Given our gene location measurements, and also our Dl grad widths and
% shapes, what amplitudes are needed to make this happen? To do this, we
% will gather all the pertinent data in this section, then create the 1x,
% 2x, and 4x Dl gradients as functions of sigma. In the next section, we
% will perform least squares analysis.
% =========================================================================
% {

%
% Load data.
%
load sig_fixed
load genes_fixed

%
% Shape of wild-type Dl and 1xdl. Note that 4xdl has the same shape as wt,
% but a different sigma.
%
b = 0.11;
m = -0.1;
Dlwt = @(x,sig)exp(-x.^2/2./sig.^2) + b + m*x;
load Dl125_avg t

s = linspace(0,1,151)';
t1 = flipud(t(1:151)); % split the 1x dl gradient in half
t2 = t(151:end);
t = 0.5*(t1 + t2); % average
t = mDU(t); % and normalize

%
% correct the 1xdl gradients for their non-gaussian shape
%
sig1x0 = mean(sig1x);
Dl1x = @(x,sig)interp1(s,t,x./(sig/sig1x0),'pchip',0); 

%}

%% ========================================================================
% Now we take our gathered data and perform the LSQ analysis.
% =========================================================================
% {

%
% Calc means and SEM's of our data
%
X = {sD_sna1x,sD_snawt,sD_sna4x
	sV_sog1x,sV_sogwt,sV_sog4x
	sD_sog1x,sD_sogwt,sD_sog4x};
SIG = {sig1x,sigwt,sig4x};
SEMx = cellfun(@std,X)./sqrt(cellfun(@length,X));
SEMsig = cellfun(@std,SIG)./sqrt(cellfun(@length,SIG));

X = cellfun(@mean,X);
SIG = cellfun(@mean,SIG);

%
% Calc Dl concentrations at the different positions and genotypes
%
DL = [Dl1x(X(:,1),SIG(1)) Dlwt(X(:,2),SIG(2)) Dlwt(X(:,3),SIG(3))];

%
% calc derivatives of Dl wrt x and to sig
%
delt = 1e-2;
dDl1xdx = (Dl1x(X(:,1)+delt,SIG(1)) - Dl1x(X(:,1),SIG(1)))/delt;
dDlwtdx = (Dlwt(X(:,2)+delt,SIG(2)) - Dlwt(X(:,2),SIG(2)))/delt;
dDl4xdx = (Dlwt(X(:,3)+delt,SIG(3)) - Dlwt(X(:,3),SIG(3)))/delt;
dDl1xdsig = (Dl1x(X(:,1),SIG(1)+delt) - Dl1x(X(:,1),SIG(1)))/delt;
dDlwtdsig = (Dlwt(X(:,2),SIG(2)+delt) - Dlwt(X(:,2),SIG(2)))/delt;
dDl4xdsig = (Dlwt(X(:,3),SIG(3)+delt) - Dlwt(X(:,3),SIG(3)))/delt;

%
% Put these together to get the variability for our LSQ
%
s = sqrt([dDl1xdx dDlwtdx dDl4xdx].^2.*SEMx.^2 + ...
	[dDl1xdsig dDlwtdsig dDl4xdsig].^2.*repmat(SEMsig.^2,3,1));

%
% Now, construct the matrix and the target
%
A1x = [sum((DL(:,1)./s(:,1)).^2) 0 -(DL(:,1)./s(:,1).^2)'];
A4x = [0 sum((DL(:,3)./s(:,3)).^2) -(DL(:,3)./s(:,3).^2)'];
Asna = [(DL(1,[1,3])./s(1,[1,3]).^2) -sum(1./s(1,:).^2) 0 0];
Asogv = [(DL(2,[1,3])./s(2,[1,3]).^2) 0 -sum(1./s(2,:).^2) 0];
Asogd = [(DL(3,[1,3])./s(3,[1,3]).^2) 0 0 -sum(1./s(3,:).^2)];

B = [0; 0; -(DL(:,2)./s(:,2).^2)];
P = [A1x; A4x; Asna; Asogv; Asogd] \ B;

%
% Results:
%
alpha1x = P(1);
alpha4x = P(2); 
theta_sna = P(3);
theta_sogv = P(4);
theta_sogd = P(5);

%
% The above results generally make sense. First, alpha1x is about 50%, but
% not quite. I think my by-eye estimate was 70%, and this gives 60%.
% Second, sogv and sna have a high threshold...and they are very close to
% each other. Third, sogd has a very low threshold.
%
% What doesn't make sense: alpha4x is less than one (90%). That's a
% bit off, but I'm ok with it, since there are undoubtedly errorbars on
% these estimates. At any rate, I will now use these estimates to create
% a graph that shows the minimal perturbation of gene expression
%
x = linspace(0,1,101)';
figure
plot(x,alpha1x*mDU(Dl1x(x,SIG(1))))
hold on
plot(x,mDU(Dlwt(x,SIG(2))))
plot(x,alpha4x*mDU(Dlwt(x,SIG(3))))
plot(xlim,theta_sna*[1 1],'k :')
plot(xlim,theta_sogv*[1 1],'g :')
plot(xlim,theta_sogd*[1 1],'k :')
set(gcf,'paperpositionmode','auto')
set(gca,'fontsize',16,'xtick',0:0.2:1,'ytick',0:0.2:1)

alpha1xstar = alpha1x;
alpha4xstar = alpha4x;

%}

%% ========================================================================
% Here we calculate the objective function for several values of amplitudes
% =========================================================================
% {

w = 1./s.^2;
D = sum(w,2);
Nhat = w.*DL;

%
% Double for-loop
%
n1x = 50;
n4x = 50;
Alpha1x = linspace(0.3,1.2,n1x)';
Alpha4x = linspace(0.8,1.6,n4x)';
alpha2x = 1;

F = zeros(n1x,n4x);
for i = 1:n1x
	alpha1x = Alpha1x(i);
	for j = 1:n4x
		alpha4x = Alpha4x(j);
		
		
		%
		% Now, calc the best-fit theta for each boundary
		%
		Alpha = [alpha1x; alpha2x; alpha4x];
		N = Nhat*Alpha;
		theta = N./D;
		
		%
		% Now that we have the best-fit theta, we can calculate our
		% objective function
		%
		f1 = ((alpha1x*DL(:,1) - theta)./s(:,1)).^2 + ...
			((alpha2x*DL(:,2) - theta)./s(:,2)).^2 + ...
			((alpha4x*DL(:,3) - theta)./s(:,3)).^2;
		f = sum(f1);
		
		F(i,j) = f;

	end
end

%
% Calc the best-fit theta for each boundary, at the best-fit alpha's
%
Alpha = [alpha1xstar; alpha2x; alpha4xstar];
N = Nhat*Alpha;
theta = N./D;

%
% Calculate fstar, which is the value of f at the best-fit alpha's.
%
f1 = ((alpha1xstar*DL(:,1) - theta)./s(:,1)).^2 + ...
	((alpha2x*DL(:,2) - theta)./s(:,2)).^2 + ...
	((alpha4xstar*DL(:,3) - theta)./s(:,3)).^2;
fstar = sum(f1);
F = F/fstar; % normalize by the global minimum value of f

figure
pcolor_contour(Alpha1x,Alpha4x,F',2:2:8,[],false)
plot(alpha1xstar,alpha4xstar,'. r')
set(gca,'fontsize',24,'xtick',0.4:0.2:1.2)

%}

function pcolor_contour(X,Y,Z,thresh,colors,vlog,Zmin,Zmax,colorbarticks)
% function pcolor_contour(X,Y,Z,thresh,colors,vlog,Zmin,Zmax,colorbarticks)
%
% X,Y: obvious
% Z: Z data used for pcolor and contour plots. If you really want to have
%	the pcolor and contour plot two different sets of Z data, let Z  be a
%	cell variable with the first element being the Z data for pcolor, and
%	the second being the Z data for contour.
% thresh (optional): row vec that specifies the values for the contours. If
%  left off, the contours will automatically be the ticks on the colorbar
% colors (optional): cell array for the colors of the contours, corresp to
%	thresh. If left off, all contours will be red.
% vlog: (optional) whether you want the x and y axes to be on the log
%	scale. Can be a 2-elt logical vector. If so, then the first elt
%	specifies whether the x axis should be in log scale, and the second for
%	whether the y axis should be in log scale. Senseless if both are true
%	(in that case, just make it a 1-elt scalar true).
% Zmin,Zmax (optional): the max and min to show for the pcolor. If values
%	are above the max, they'll show as yellow (highest caxis color for
%	'parula' colormap), and if less than min, they'll show as dark blue. If
%	not specified, then caxis limits will be default.
% colorbarticks (optional): the tick mark labels for the colorbar. Should
%	be passed as a row vector. If passed as an empty variable, no colorbar
%	will be shown.

%
% Check to see if the Z data for the pcolor are different from those for
% the contour
%
if iscell(Z)
	Zp = Z{1};
	Zc = Z{2};
else
	Zp = Z;
	Zc = Z;
end

%
% Make pcolor plot
%
% figure('paperpositionmode','auto')
pcolor(X,Y,Zp)
shading interp
hold on

%
% Check caxis limits, if asked for
%
CAXIS = caxis;
if exist('Zmin','var') && ~isempty(Zmin) && ~isnan(Zmin) && isnumeric(Zmin)
	if exist('Zmax','var') && ~isnan(Zmax) && isnumeric(Zmax) && ~isempty(Zmax)
		caxis([Zmin Zmax])
	else
		caxis([Zmin CAXIS(2)])
	end
elseif exist('Zmax','var') && ~isempty(Zmax) && ~isnan(Zmax) && isnumeric(Zmax)
	caxis([CAXIS(1) Zmax])	
end

%
% Set colorbarticks, if asked for
%
if exist('colorbarticks','var')
	if ~isempty(colorbarticks)
		colorbar
		h = get(gcf,'children');
		for i = 1:length(h)
			if strcmp(h(i).Tag,'Colorbar')
				set(h(i),'YTick',colorbarticks)
			end
		end
	end
else
	colorbar
end

%
% Change axis (or axes) to log scale, if asked for.
%
if exist('vlog','var') && length(vlog) == 1 && vlog
	set(gca,'xscale','log','yscale','log')
elseif exist('vlog','var') && length(vlog) == 2
	if vlog(1)
		set(gca,'xscale','log')
	elseif vlog(2)
		set(gca,'yscale','log')
	end
end

%
% Threshold check
%
if isempty(thresh)
	h = get(gcf,'children');
	for i = 1:length(h)
		if strcmp(h(i).Tag,'Colorbar')
			thresh = get(h(i),'YTick');
			break
		end
	end
	thresh1 = thresh; 
elseif length(thresh) == 1
	thresh1 = [thresh thresh];
elseif size(thresh,1) > 1 && size(thresh,2) > 1
	error('thresh param has too many dimensions. See if you passed two Z''s')
else
	thresh1 = thresh; 
end

%
% Colors check
%
if isempty(colors)
	colors = {'r'};
elseif ~iscell(colors) && length(colors) == 1
	colors = {colors}; 
end
if length(colors) == 1
	colors = repmat(colors,1,length(thresh)); 
end
leg = cell(size(thresh));
Vthresh = false(size(thresh));

%
% Plot contours
%
c = contourc(X,Y,Zc,thresh1);
if ~isempty(c)
	[X1,Y1,C] = extractcontour(c);
	for k = 1:length(Y1)
		vthresh = C(k) == thresh;
		plothand = plot(X1{k},Y1{k},colors{vthresh});
		if isempty(leg{vthresh})
			leg{vthresh} = num2str(thresh(vthresh));
			plothandles(vthresh) = plothand; %#ok<AGROW>
			Vthresh = Vthresh | vthresh;
		end
	end
	legend(plothandles(Vthresh),leg(Vthresh))
end


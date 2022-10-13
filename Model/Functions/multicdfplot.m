function varargout = multicdfplot(y,varargin)
%plots cumulative distribution function of "x"
%
%function varargout = cdfplotDU(x,varargin)
%
% This function takes your variable "x" and performs the following:
%	"plot(sort(x),(0:(length(x)-1))/(length(x)-1))"
%
%	NOTE: the default linewidth is 2.
%
% "x": is a col vec or array of values.  Arrays will be assumed to have the
%	relevant data sets going down their columns.
% 
% Optional argument varargin can consist of any set of additional arguments
% that you could normally pass to Matlab's "plot" function.  A special case
% is when the first element of varargin is a scalar integer, in which case
% this function will use that for the figure number.  The remaining
% arguments then in varargin will be passed as options to the "plot"
% function.
%
% varargout will be {h}, the handle(s) of the plot, if asked for.

%
% Unpacking varargin.
%
nArg = size(varargin,2); 
iArg = 1; 

if nArg >= iArg && (isnumeric(varargin{iArg}) || isgraphics(varargin{iArg}))
	fignum = varargin{iArg};
	if nArg > iArg
		varargin = varargin(2:end);
	else
		varargin = {};
	end
end%, iArg = iArg + 1;

if exist('fignum','var')
	figure(fignum)
else
	figure
end

len = length(y);

for i=1:len
    x = y(:,i);
    
    if any(isnan(x))
        disp('Warning: NaN values ignored')
    end
    x(isnan(x)) = [];
    if length(x) ~= 1
        h(i) = plot(sort(x),(0:(size(x,1)-1))/(size(x,1)-1),varargin{:});
    else
        h(i) = plot([x x],[0 1],varargin{:});
    end
end

for i = 1:length(h)
	set(h(i),'Linewidth',2)
end
if nargout > 0
	varargout = {h};
end
end

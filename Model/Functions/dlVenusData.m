function [C,dC] = dlVenusData(options,varargin)

% dlVenusData(options) - a function that creates a column vector of the
% dl-Venus data from Reeves et al., 2012, useful for calculating error
%
% C : the column vector of data
% dC: the weights of each data point
%
% options: 'nuclear' - nuclear dl fluorescence
%          'totalDl' - total dl fluorescence

load('Mat/hybrid_embryo.mat')
switch options
    case 'nuclear'
        A = data.A;
        B = data.B;
        dA = data.dA;
        dB = data.dB;
    case {'totalDl','totalDorsal','total'}
        A = data.A2;
        B = data.B2;
        dA = data.dA;
        dB = data.dB;
end

args = struct('M',[19 26 36 51],'yesplot',false);
if nargin > 1
    if mod(length(varargin),2)==0
        for i = 1:2:length(varargin)
            args.(varargin{i})=varargin{i+1};
        end
    else
        error('Varargin requires name-value pairs. Some options not specified')
    end
end

%% NC 11
% N11T = t(1:16);
N11A = A(1:16);
N11B = B(1:16);
dA11 = dA(1:16);
dB11 = dB(1:16);

sig = 0.15;
m = -0.07;
M = args.M(1);
x = linspace(0,1,M);
N11D = repmat(N11A,1,M).*exp(-repmat(x,length(N11A),1).^2/2/sig^2) + ...
    repmat(N11B,1,M) + m*N11A*x;
N11D = fliplr(rot90(N11D,3));

dC11 = exp(-repmat(x,length(N11A),1).^2/2/sig^2).*repmat(dA11,1,length(x)) + ...
	m*dA11*x + m*repmat(x,length(dA11),1) + repmat(dB11,1,length(x));

%% NC 12
% N12T = t(34:59);
N12A = A(34:59);
N12B = B(34:59);
dA12 = dA(34:59);
dB12 = dB(34:59);
M = args.M(2);
x = linspace(0,1,M);
N12D = repmat(N12A,1,M).*exp(-repmat(x,length(N12A),1).^2/2/sig^2) + ...
    repmat(N12B,1,M) + m*N12A*x;
N12D = fliplr(rot90(N12D,3));
dC12 = exp(-repmat(x,length(N12A),1).^2/2/sig^2).*repmat(dA12,1,length(x)) + ...
	m*dA12*x + m*repmat(x,length(dA12),1) + repmat(dB12,1,length(x));

%% NC 13
% N13T = t(77:126);
N13A = A(77:126);
N13B = B(77:126);
dA13 = dA(77:126);
dB13 = dB(77:126);
M = args.M(3);
x = linspace(0,1,M);
N13D = repmat(N13A,1,M).*exp(-repmat(x,length(N13A),1).^2/2/sig^2) + ...
    repmat(N13B,1,M) + m*N13A*x;
N13D = fliplr(rot90(N13D,3));
dC13 = exp(-repmat(x,length(N13A),1).^2/2/sig^2).*repmat(dA13,1,length(x)) + ...
	m*dA13*x + m*repmat(x,length(dA13),1) + repmat(dB13,1,length(x));

%% NC 14
% N14T = t(149:332);
N14A = A(149:332);
N14B = B(149:332);
dA14 = dA(149:332);
dB14 = dB(149:332);

M = args.M(4);
x = linspace(0,1,M);
N14D = repmat(N14A,1,M).*exp(-repmat(x,length(N14A),1).^2/2/sig^2) + ...
    repmat(N14B,1,M) + m*N14A*x;
N14D = fliplr(rot90(N14D,3));
dC14 = exp(-repmat(x,length(N14A),1).^2/2/sig^2).*repmat(dA14,1,length(x)) + ...
	m*dA14*x + m*repmat(x,length(dA14),1) + repmat(dB14,1,length(x));

%% Fluorescence data as a column vector
C = [N11D(:);N12D(:);N13D(:);N14D(:)];
dC = [dC11(:);dC12(:);dC13(:);dC14(:)];
% totalDorsal_int = [];





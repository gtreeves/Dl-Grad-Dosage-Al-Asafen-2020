% script_run_liveanalysis_GFP
%
% This script analyzes Hadel's Dorsal-mGFP 1X, 2X and 4X embryos.

clear
close all

pthbase = 'Datadryad_live'; % must create the correct base to the 
% path to the live imaging folders on the local computer.

pth = {'2018-01-26'
	'2018-01-30'
	'2018-02-01'
	'2018-02-02'
	'2018-02-07' 
	'2018-02-11'
	'2018-02-14'
	'2018-02-24'
	'2018-05-24 Live'
	'2018-06-01'
	'2018-06-02'
	};

%
% Image parameters
%	
ring_width = 18.36; % microns
nuc_ch = 1; % nuclear channel. EDIT THIS FOR EVERY DIFFERENT FILE RUN
np_ch = 2; % nuclear protein channel
nt = 60;
z0 = 150; % how deep the slice is
h = 0.25;

%
% Loop over every path
%
Soln = [];
for iii = 1:length(pth)
	pth2 = [pthbase,pth{iii}];
	[filenames,filenamesshort] = readdir2(pth2,'lsm');	
	
	%
	% Find the calibration lsm file, store its name as a separate entity, and
	% remove it from the list of files we will do analysis on.
	%
	v = strfindiDU(filenamesshort,'calibration.lsm');
	calibration_filename = filenames{v};
	filenames(v) = [];
    nLSM = length(filenames);
	
	%
	% Calculate the calibration
	%
	[IM,lsminf1,lsminf2] = dosage_imread(calibration_filename);
	calib_mean = mean(mean(double(IM(:,:,2))));
	calib_LP = lsminf2.ScanInfo.POWER;
	
	%
	% Run analysis on all lsm files (except the calibration)
	%
	for ii = 1:nLSM
		filename = filenames{ii};
		
		if any(strfind(filename,'1 copy')) || any(strfind(filename,'1X')) 
			genotype = 'dl1x';
		elseif  any(strfind(filename,'2 copy')) || any(strfind(filename,'2X'))
			genotype = 'dl2x'; 
		elseif  any(strfind(filename,'4 copy')) || any(strfind(filename,'4X')) 
			genotype = 'dl4x';
		else
			genotype = 'unknown';
		end
		
		%
		% Reading in the lsm data
		%
		[IM,lsminf1,lsminf2] = dosage_imread(filename);
		if strcmp(filename,[pthbase,'2018-02-01',filesep,...
				'2018-02-01 2 copy dorsal embryo 01.lsm'])
			IM(:,:,:,201:end) = []; % This file ran too long.
		elseif strcmp(filename,[pthbase,'2018-01-26',filesep,...
				'2018-01-26 4 copy embryo 01.lsm'])
			IM(:,:,:,62:end) = []; % This file also ran too long.
		elseif strcmp(filename,[pthbase,'2018-05-25',filesep,...
				'2018-05-25 2 copy dorsal embryo 17.lsm'])
			IM(:,:,:,end) = []; % Final frame of that embryo was incomplete
		elseif strcmp(filename,[pthbase,'2018-06-01',filesep,...
				'4 copy dorsal embryo 26.lsm'])
			IM(:,:,:,69:end) = []; % Final frames of that embryo was
			% incomplete
		end
		
		%
		% Metadata
		%
		scalings = 1e6*[lsminf2.VoxelSizeX,lsminf2.VoxelSizeY,...
			lsminf2.VoxelSizeZ]; % microns/pixel
		Yhatmax = round(ring_width/scalings(1)); % This is the width of the ring
		% (radius_outer - radius_inner) in pixels.  Default width: 18.36 microns.
		dt = lsminf2.TimeStamps.AvgStep;
		rho = round(scalings(3)/scalings(1)); % hope we don't have to round much.
		LP = [lsminf2.ScanInfo.POWER{1} lsminf2.ScanInfo.POWER{2}];
		[H,W,numch,D] = size(IM);
		
		%
		% Subtracting background.  We assume the background intensity level (i.e.,
		% the true black level) for each channel is equal to the mode of
		% intensities seen in a slice.
		%
		bg = zeros(numch,1);
		bgsig = zeros(numch,1);
		for i = 1:numch
			slice = double(IM(:,:,i,1));
			[n,x] = hist(slice(:),0:65535);
			[A,k] = max(n);
			bg(i) = x(k);
			X = x(max(k-4,1):min(k+4,65536));
			Y = n(max(k-4,1):min(k+4,65536));
			SIG = (X-bg(i))./sqrt(-2*log(Y/A));
			bgsig(i) = sqrt(meanDU(SIG'.^2));
			IM(:,:,i,:) = imsubtract(IM(:,:,i,:),bg(i));
		end
				
		%
		% First we find the periphery of the embryo and detect the intensity of
		% each channel in quadrilaterals as you go around the periphery.  Next, we
		% perform calculations on the nuclei/nuclear proteins/nuclear dots, if
		% these things are present.  To do these things, we go through a loop that
		% takes each slice individually.
		%		
		arc = zeros(D,1);
		Xp = cell(D,1); Yp = cell(D,1);
		X1 = cell(D,1); Y1 = X1; S = X1;
		Nuc = X1; Std_nuc = X1;
		Nuc_protein = X1; Std_nuc_protein = X1;
		Intron = X1; Std_intron = X1;
		U1 = X1;
		for i = 1:D
			
			I = IM(:,:,:,i);
			
			%
			% border and perimiter.
			%
			yesplot = true;
			[xp,yp] = borderFinder(I,h,yesplot,nt,[],[],filename,i);
			Xp{i} = xp; Yp{i} = yp;
			arc(i) = sum(sqrt(diff(xp).^2 + diff(yp).^2));
			
			
			%
			% Now we detect the nuclei in each slice, given a nuclear channel.
			%
			if ~isempty(nuc_ch)
				%
				% The function "find_nuclei" detects the nuclei (see the function
				% description for more info on how) and returns the pixel locations
				% of the nuclei, the locations of the centroids of the nuclei,
				% these locations in terms of pseudo-arclength, the perimeter of
				% the embryo, and a nuclear mask image.
				%
				stage = '';
				[nucstats,xnuc,ynuc,snuc,w,mask,u1] = ...
					find_nuclei(I(:,:,nuc_ch),...
					xp,yp,scalings,Yhatmax,yesplot,stage,filename,i);
				% 		U1{i} = u1;
				U1{i} = mask;
				X1{i} = xnuc;
				Y1{i} = ynuc;
				S{i} = 2*snuc/w - 1;
				
				%
				% Obtaining the intensity of the nuclear protein channels, as well
				% as the nuclear channel itself, for each nucleus.
				%
				[np,stdnp] = nuclearintensity(I(:,:,[nuc_ch np_ch]),nucstats);
				Nuc{i} = np(:,1);
				Std_nuc{i} = stdnp(:,1);
				if ~isempty(np_ch)
					Nuc_protein{i} = np(:,2:end);
					Std_nuc_protein{i} = stdnp(:,2:end);
				elseif i == D
					Nuc_protein = NaN; Std_nuc_protein = NaN;
				end
				
			else
				X1 = NaN; Y1 = NaN; S = NaN;
				Nuc = NaN; Std_nuc = NaN;
				Nuc_protein = NaN; Std_nuc_protein = NaN;
				Intron = NaN; Std_intron = NaN;
			end
		end
		
		%
		% Transforming data on the nuclear protein into the scaled data
		%
		nuc = cell2mat(Nuc);
		std_nuc = cell2mat(Std_nuc);
		R = cell(D,1); Std_R = cell(D,1);
		for i = 1:D
			%
			% An explanation of the following calculation.
			%
			Nuc_norm = Nuc{i}./median(Nuc{i}); % Dividing by this corrects
			% for changes in intensity around the circumference of the
			% embryo, at frame i. It is the nuclear intensity, which is
			% supposed to be uniform, normalized by the median nuclear
			% intensity (so that dividing by this vector does not bring our
			% nuc protein intensity down to order 1.
			Timecourse_norm = median(nuc)/LP(1)/100; % Dividing by this
			% scalar corrects for an overall brighter timecourse in the
			% nuclei. This could be an issue if, say, one embryo was taken
			% at z = 160 and another at z = 140, by accident. The overall
			% nuclear intensity should be the same across every embryo
			% imaged.
			Calib_norm = calib_mean/1.75e4; % Dividing by this scalar 
			% corrects for day-to-day laser power changes. The calib_mean
			% is the mean intensity of the calibration image taken on a
			% given day.
			
			R{i} = Nuc_protein{i}/LP(2) ./ Nuc_norm / Timecourse_norm / Calib_norm;
			
			Std_R{i} = R{i}.*sqrt((Std_nuc_protein{i}./Nuc_protein{i}).^2 + ...
				repmat((Std_nuc{i}./Nuc{i}).^2,1,size(Nuc_protein{i},2)));
			
		end
		
		%
		% "data" remarks:
		%
		%%
		data.filename = filename;
		data.genotype = genotype;
		data.metadata.lsminf1 = lsminf1;
		data.metadata.lsminf2 = lsminf2;
		data.metadata.scalings = scalings;
		data.metadata.rho = rho;
		data.metadata.Yhatmax = Yhatmax;
		data.metadata.nt = nt;
		data.metadata.bg = bg;
		data.metadata.std_bg = bgsig;
		data.H = H;
		data.W = W;
		data.D = D;
		data.dt = dt;
		data.z0 = z0;
		
		data.Lbar = mean(arc)/2*scalings(1);
		data.s_mid = NaN;
		data.S = S;
		data.R = R;
		data.Std_R = Std_R;
		
		data.metadata.LP = LP;
		data.metadata.calib_mean = calib_mean;
		data.metadata.calib_LP = calib_LP;
		data.metadata.w = mean(arc);
		data.metadata.arc = arc;
		data.metadata.L = arc/2*scalings(1);
		data.metadata.Xp = Xp;
		data.metadata.Yp = Yp;
		data.metadata.X = X1;
		data.metadata.Y = Y1;
		data.metadata.S = S;
		data.metadata.Nuc = Nuc;
		data.metadata.Std_nuc = Std_nuc;
		data.metadata.Nuc_protein = Nuc_protein;
		data.metadata.Std_nuc_protein = Std_nuc_protein;
		
		
		data.nucprotein_names = {};
		data.A = NaN;
		data.B = NaN;
		data.M = NaN;
		data.Sig = NaN;
		data.dA = NaN;
		data.dB = NaN;
		data.dM = NaN;
		data.dSig = NaN;
		data.gof = NaN;
		
		save([data.filename(1:end-4),'_data.mat'],'data');
		
		disp(['ii = ',num2str(ii),' out of ',num2str(nLSM),...
			';   iii = ',num2str(iii),' out of ',num2str(length(pth))])
	end
end


clk1 = clock;
clk = [num2str(clk1(1)),'-',num2strDU(clk1(2),2),'-',num2strDU(clk1(3),2),'_',...
	num2strDU(clk1(4),2),'-',num2strDU(clk1(5),2),'-',num2strDU(round(clk1(6)),2)];

save(['Mat/',clk,'_Soln'],'Soln')

%}

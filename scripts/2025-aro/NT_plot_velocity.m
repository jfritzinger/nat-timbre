%% NT_velocity

%% calculate_velocity_stimulus
clear

%% Load in stimulus

if ismac
	fpath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data';
else
	fpath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\';
end
tuning = readtable(fullfile(fpath, 'Tuning.xlsx')); % Load in tuning

target = 'Bassoon';
listing = dir(fullfile(fpath, 'waveforms', '*.wav'));
target_WAV = arrayfun(@(n) contains(listing(n).name, target), 1:numel(listing), 'UniformOutput', false);
wav_nums =  find(cell2mat(target_WAV));

d = dir(fullfile(fpath,'waveforms', '*.wav'));
all_files = sort({d.name});
nfiles = length(wav_nums);
wav_npts = zeros(1,nfiles);
wav_data = cell(1,nfiles);

for i = 1:nfiles
	files{1,i} = all_files{wav_nums(i)};
end

% Sort by frequency of pitch
index = [];
note_names = extractBetween(files, 'ff.','.');

for ii = 1:nfiles % Find index of each note in tuning spreadsheet
	index(ii) = find(strcmp(note_names(ii), tuning.Note));
end
pitch_order = tuning.Frequency(index); % Get freqs of each note
[~, order] = sort(pitch_order); % Sort freqs
F0s = pitch_order;
files = files(order);
F0s = F0s(order);

itype = 1;
criteria2 = [2 2 2 2 2 2 2 2 2 2 ...
	2 2 2 2 1.5 2 2 2 3 2 ...
	1.6 2 2 2 2 2 2 2 2 3 ...
	2 2 2 2 2 2 2 2 2 2];
for ind =  21 %20 %1:nfiles
	% Originally 21, 20
	% Good options: 11

	target_F0 = F0s(ind);
	target_file = fullfile(fpath,'waveforms', files{ind});

	%% Paul's spectrogram

	% Creates filename for analyzed data
	if ismac
		savepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/decomped';
	else
		savepath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\decomped';
	end
	savefile = sprintf('%s_F0_%0.0f_decomp.mat', target, target_F0);

	temp = load(fullfile(savepath, savefile));
	decomp_info = temp.decomp_info;

	% Sets time vector & max frequency
	[data,fs] = audioread(target_file);
	dur = length(data)/fs; %in sec
	t = (1/fs):(1/fs):dur;
	t = t(1:length(data)); % not sure why t and data are different lengths, this remedies that

	%%%%%% Could change this to be a specific # of harmonics each time
	max_freq = target_F0*25; %3500; % Determines maximum frequency to analyze
	if max_freq>fs/2
		max_freq = fs/2-1;
	end
	%%%%%%

	% get decomp variables out of struct
	starti = decomp_info(1).start_inds{1};
	endi = decomp_info(1).end_inds{1};
	num_parts = decomp_info(1).num_parts;

	% construct "not-spectrogram" spectrogram using phase first, append
	% waveforms of cosines with appropriate phase into the same array
	phase_spect = zeros(length(decomp_info.freq{1}),length(t));
	for icomp = 1:length(decomp_info.freq{1})
		A2 = decomp_info.A{1}(icomp);
		f2 = decomp_info.freq{1}(icomp);
		phi2 = decomp_info.phi{1}(icomp);
		phase_spect(icomp,:) = (A2*cos(2*pi*t*f2 + phi2));
	end

	% manually unwrap phase so that it is monotonically decreasing
	for parti = 1:num_parts
		if parti == 1
			data_part = data(starti(parti):endi(parti));
			A = decomp_info(parti).A{:};
			freq = decomp_info(parti).freq{:};
			phi = decomp_info(parti).phi{:};
			uwphi = [];
			for iphi = 1:length(phi)
				if iphi == 1
					uwphi(iphi) = phi(iphi);
				else
					temp = phi(iphi);
					while temp > uwphi(iphi-1)
						temp = temp - 2*pi;
					end
					uwphi(iphi) = temp;
				end
			end
		end

		% Reconstruct signal using these phases
		temp_pin = zeros(1,length(t(starti(parti):endi(parti))));
		for iharm = 1:length(freq)
			% without magnitude
			temp_pin = temp_pin + cos(2*pi*t(starti(parti):endi(parti))*freq(iharm) + phi(iharm));

			% with magnitude applied
			% temp_pin = temp_pin + A(iharm)*cos(2*pi*t(starti(parti):endi(parti))*freq(iharm) + phi(iharm));
		end
		pin = temp_pin;

		% Generate spectrogram
			if itype == 1
				window = round(1/decomp_info(parti).F0_actual*60*1000); %round(decomp_info(parti).F0_actual*6); %600;
				ov = window-10; %round(window*0.9833); %590; % or possibly try minus 10 instead?
				[sg,Ftmp,Ttmp] = spectrogram(pin,window,ov,[],fs,'yaxis');
			elseif itype == 2
				window = round(1/decomp_info(parti).F0_actual*fs/2);
				ov = round(0.95*window);
				[sg,Ftmp,Ttmp] = spectrogram(pin,window,ov,[],fs,'yaxis');
			else
				window = round(1/decomp_info(parti).F0_actual*fs/2);
				ov = round(0.95*window);
				[sg,Ftmp,Ttmp] = spectrogram(data,window,ov,[],fs,'yaxis');
			end

		Ttmp_part(parti) = {Ttmp};
		Ftmp_part(parti) = {Ftmp};
		spect(parti) = {20*log10(abs(sg))};
		fi = 1;
		while Ftmp(fi) < max_freq
			fi = fi+1;
		end

	end

	for parti = 1:num_parts
		this_spect = spect{parti};
		Ttmp = Ttmp_part{parti};
		Ftmp = Ftmp_part{parti};
	end

	%% Finding peaks using a threshold

	F0 = F0s(ind);
	period_lim = (1/F0)*6;
	ylimits = [0 max_freq]/1000;
	xlimits = [0.02 period_lim+0.02];

	% Cut down to avoid onsets/offsets
	Fs = fs;
	this_spect = this_spect(1:fi,:);
	Ftmp = Ftmp(1:fi);
	t = Ttmp;

	% Using findpeaks function
	num_harms = floor(max_freq/F0);
	harms = F0:F0:max_freq;
	peaks_all = NaN(size(num_harms,1), size(this_spect,2));
	spec_fs = length(Ttmp)/Ttmp(end);
	npts = floor(1/F0*spec_fs);
	for j1 = 1:num_harms

		[~, k1] = min(abs(Ftmp-harms(j1)));
		[pks,loc] = findpeaks(this_spect(k1,:),'MinPeakDistance',npts/1.5);
		P{k1} = [pks; loc];
		peaks = NaN(1,length(t));
		peaks(loc) = Ftmp(k1); %CF;
		peaks_all(j1, :) = peaks;
	end

	% Get an array of timing and frequency information
	peaks_all(isnan(peaks_all)) = 0;
	[rows, cols, values] = find(peaks_all);
	timing = t(cols);
	freqs = values;
	freqs2 = rows;


	% Plot
	% t = Ttmp;
	% ax1 = nexttile([1, 2]);
	% hold on
	% scatter(t, peaks_all./1000,30, 'filled' , 'MarkerEdgeColor','k',...
	% 	'MarkerFaceColor','k')
	% xlabel('Time (ms)')
	% title('Peaks')
	% set(gca, 'FontSize', 16)
	% ylim(ylimits)
	% xlim(xlimits)
	% grid on

	% 'Unwrap' properly
	%figure
	%hold on
	criteria = (1/F0)/criteria2(ind);
	time_example = zeros(1, max(freqs2));
	for iCF = 1:max(freqs2)
		indx = find(iCF==freqs2);
		time = timing(indx(8));

		if iCF ~= 1
			max_time = time_last+criteria;
			min_time = time_last-criteria;
			if time > max_time
				time = timing(indx(7));
			elseif time < min_time
				time = timing(indx(9));
			end
		end
		%scatter(time,harms(iCF), 'filled');
		time_last = time;
		time_example(iCF) = time;
	end

	fontsize = 30;

	% Shift so all unwrap properly
	max_f = max(freqs2);
	start_t = min(time_example);
	t_avg = time_example-start_t;
	% timing = timing - start_t;

	% Calculate df/dt for each nearest neighbor
	dt = diff(t_avg*1000);
	df = diff(harms(1:max_f)/1000);
	dt_zero_ind = find(dt==0);
	v = df./dt;
	v(dt_zero_ind) = 0;

	%% Full neurogram
	figure('Position',[100,100,1400,500])
	tiledlayout(1, 7, 'TileSpacing','compact')
	ax(1) = subplot(1, 4, 1);
	s = surf(Ttmp,Ftmp./1000,this_spect);
	s.EdgeColor = 'none';
	%shading interp
	hold on
	view(0,90)
	xlim(xlimits)
	title('Spectrogram')
	xlabel('Time (s)')
	ylabel('CFs (kHz)')
	ylim(ylimits)
	set(gca, 'FontSize', fontsize)
	if itype == 1
		ax.CLim = [35, max(this_spect(1:fi,:), [], 'all')];
		colormap(gray)
	elseif itype == 2
		ax.CLim = [25, max(this_spect(1:fi,:), [], 'all')-3];
	else
		ax.CLim = [-35, max(this_spect(1:fi,:), [], 'all')];
	end
	%colormap(gray)
	plot3(t_avg+0.026, harms(1:max_f)/1000, repmat(60, length(t_avg), 1), 'r', 'LineWidth',4)
	%plot3(t_avg+0.028, harms(1:max_f)/1000, repmat(60, length(t_avg), 1), 'r', 'LineWidth',4)

	ax(2) = subplot(1, 4, 2);
	scatter(v,(harms(1:max_f-1)+F0/2)/1000, 'filled')
	hold on
	plot(v,(harms(1:max_f-1)+F0/2)/1000)
	hold on
	xlabel('Velocity (kHz/ms)')
	title('Velocity')
	set(gca, 'FontSize', fontsize)
	xline(0, 'k')
	max_vel = max(abs(v));
	xlim([-1*max_vel-0.1 max_vel+0.1])
	grid on
	v_harms = harms(1:max_f-1)+F0/2;
	ylim(ylimits)

end

%% Add RVF 


if ismac
	fpath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data';
else
	fpath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\';
end
tuning = readtable(fullfile(fpath, 'Tuning.xlsx')); % Load in tuning
target = 'Bassoon';

% Get all .wav files containing the target instrument name
listing = dir(fullfile(fpath, 'waveforms', ['*' target '*.wav']));
files = {listing.name};

% Extract note names and find corresponding frequencies
note_names = extractBetween(files, 'ff.', '.');
[~, index] = ismember(note_names, tuning.Note);
F0s = tuning.Frequency(index);

% Sort files and frequencies by pitch
[F0s, order] = sort(F0s);
files = files(order);

% Initialize variables for later use (if needed)
nfiles = numel(files);
wav_npts = zeros(1, nfiles);
wav_data = cell(1, nfiles);

[base, datapath, ~, ppi] = getPaths();
sheetpath = 'data/2025-manuscript/data-cleaning';
spreadsheet_name = 'PutativeTable2.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);

% Example session
putative = 'R29_TT3_P5_N02';
filename = sprintf('%s.mat', putative);
load(fullfile(datapath,'neural_data', filename)), 'data';
index = find(cellfun(@(s) strcmp(putative, s), sessions.Putative_Units));
CF = sessions.CF(index);
MTF_shape = sessions.MTF{index};

% Plot NT
params_NT = data(14, 2); % 13 is oboe, 14 is bassoon
data_NT = cell(2, 1);
if ~isempty(params_NT{1})
	data_NT = analyzeNT(params_NT{1});
end

% Plot RVF
params_RVF = data{5, 1};
ax(3) = subplot(1, 4, 3);
[fig, data_RVF] = plotPhysRVF([], params_RVF, []);
set(gca, 'FontSize', fontsize)


%% Add prediction

if ismac
	savepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/decomped';
else
	savepath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\decomped';
end

vel_CF = NaN(1, nfiles);
for ii = 1:nfiles

	% Get single stimulus
	target = extractBefore(files{ii}, '.');
	target_F0 = extractBetween(files{ii}, 'ff.', '.wav');
	% savefile = sprintf('%s_F0_%s_Velocity.mat', target, target_F0{1});
	% load(fullfile(savepath, savefile),'vel_all', 'CFs_all')

	savefile = sprintf('%s_F0_%d_Velocity.mat', target, round(F0s(ii)));
	load(fullfile(savepath, savefile), 'v', 'v_harms')
	vel_all = v;
	CFs_all = v_harms;
	
	% Find 
	[~, CF_ind] = min(abs(CF-CFs_all));
	vel_CF(ii) = vel_all(CF_ind);
end

% Rearrange
[sorted_vel, ind] = sort(data_RVF.velocities);
sorted_rate = data_RVF.rate(ind);

% Interpolate (linear)
RVF_vel = -9:0.1:9;
RVF_rates = interp1(sorted_vel, sorted_rate, RVF_vel, 'linear');

% Match each velocity to a rate 
for ii = 1:nfiles

	% Find rate most similar to velocity
	this_vel = vel_CF(ii);
	[~, index] = min(abs(this_vel-RVF_vel));
	if index == 1 
		rates_vel(ii) = RVF_rates(index+1);
	elseif index == length(RVF_vel)
		rates_vel(ii) = RVF_rates(index-1);
	else
		rates_vel(ii) = RVF_rates(index);
	end
end

% Plot natural timbre
ax(4) = subplot(1, 4, 4); 
hold on
errorbar(F0s,data_NT.rate, data_NT.rate_std/sqrt(params_NT{1}.nrep), 'LineWidth',2)
xlabel('Pitch (Hz)')
plot(F0s, rates_vel);
r = corrcoef(rates_vel, data_NT.rate);
R = r(1,2);
title('Example')
ylim([0 120])
set(gca, 'fontsize', fontsize)
set(gca, 'xscale', 'log')
xticks([60 110 220 440])
legend('Data', 'Prediction')
ylabel('Avg. Rate (sp/s)')
xlim([58 588])

%% Arrange 

set(ax(1), 'position', [0.05 0.2 0.2 0.7])
set(ax(2), 'position', [0.28 0.2 0.1 0.7])
set(ax(3), 'position', [0.47 0.2 0.23 0.7])
set(ax(4), 'position', [0.78 0.2 0.2 0.7])

%% Export 

savepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/figures/2025-aro';
set(gcf, 'Renderer', 'painters')
print('-dsvg', '-vector', fullfile(savepath,'NT_plot_velocity.svg'))

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

criteria2 = [2 2 2 2 2 2 2 2 2 2 ...
	2 2 2 2 1.5 2 2 2 3 2 ...
	1.6 2 2 2 2 2 2 2 2 3 ...
	2 2 2 2 2 2 2 2 2 2];
for ind = 1:nfiles
	% Errored: 7
	% Good: 11, 15, 18

	figure('Position',[828,20,853,350])
	tiledlayout(1, 6, 'TileSpacing','compact', 'Padding','compact')

	for itype = 1
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

				%pwelch(data_part,[],[],[],fs);
				%close
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
				altplot(iharm,:) = cos(2*pi*t*freq(iharm) + phi(iharm));
			end
			pin = temp_pin;

			% Generate spectrogram
			if itype == 1
				window = round(1/decomp_info(parti).F0_actual*60*1000); %round(decomp_info(parti).F0_actual*6); %600;
				ov = window-10; %round(window*0.9833); %590; % or possibly try minus 10 instead?
				[sg,Ftmp,Ttmp] = spectrogram(pin,window,ov,[],fs,'yaxis');
				
			elseif itype == 2
				window = round(1/decomp_info(parti).F0_actual*fs/2);
				ov = round(0.9*window);
				[sg,Ftmp,Ttmp] = spectrogram(pin,window,ov,[],fs,'yaxis');
			else
				window = round(1/decomp_info(parti).F0_actual*fs/2);
				ov = round(0.9*window);
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
		period_lim = (1/F0)*2;
		ylimits = [0 max_freq]/1000;
		xlimits = [0.02 period_lim+0.02];

		% Cut down to avoid onsets/offsets
		Fs = fs;
		this_spect = this_spect(1:fi,:);
		Ftmp = Ftmp(1:fi);
		t = Ttmp;

		% Full neurogram
		ax = nexttile([1,2]);
		% s = surf(Ttmp,Ftmp./1000,this_spect);
		% s.EdgeColor = 'none';
		% view(0,90)
		% xlim(xlimits)
		% title(sprintf('Spectrogram, F0=%0.0f', F0))
		% xlabel('Time (s)')
		% ylabel('CFs (kHz)')
		% ylim(ylimits)
		% set(gca, 'FontSize', 16)
		% %ax.CLim = [min(this_spect(1:fi,:), [], 'all'), max(this_spect(1:fi,:), [], 'all')];
		% if itype == 1
		% 	ax.CLim = [30, max(this_spect(1:fi,:), [], 'all')];
		% 	colormap(gray)
		% elseif itype == 2
		% 	ax.CLim = [15, max(this_spect(1:fi,:), [], 'all')];
		% else
		% 	%ax.CLim = [0, max(this_spect(1:fi,:), [], 'all')];
		% end

		%figure
		f_alt = F0:F0:max_freq;
		t_alt = (1/fs):(1/fs):dur;
		s = surf(t_alt,f_alt./1000,altplot);
		s.EdgeColor = 'none';
		view(0,90)
		xlim(xlimits)
		title(sprintf('Spectrogram, F0=%0.0f', F0))
		xlabel('Time (s)')
		ylabel('CFs (kHz)')
		ylim(ylimits)
		set(gca, 'FontSize', 16)
		shading interp

		% Use new
		this_spect = altplot;
		Ftmp = f_alt;
		Ttmp = t_alt;
		t = Ttmp;

		% Using findpeaks function
		num_harms = floor(max_freq/F0);
		harms = F0:F0:max_freq;
		peaks_all = NaN(size(num_harms,1), size(this_spect,2));
		spec_fs = length(Ttmp)/Ttmp(end);
		npts = floor(1/F0*spec_fs);
		for j1 = 1:num_harms
			[~, k1] = min(abs(Ftmp-harms(j1)));
			[pks,loc] = findpeaks(this_spect(k1,:)); %,'MinPeakDistance',npts); %/1.5);
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
		t = Ttmp;
		ax1 = nexttile([1, 2]);
		hold on
		scatter(t, peaks_all./1000,30, 'filled' , 'MarkerEdgeColor','k',...
			'MarkerFaceColor','k')
		xlabel('Time (ms)')
		title('Peaks')
		set(gca, 'FontSize', 16)
		ylim(ylimits)
		xlim(xlimits)
		grid on

		% 'Unwrap' properly
		figure
		hold on
		criteria = (1/F0)/6;
		time_example = zeros(1, max(freqs2));
		for iCF = 1:max(freqs2)
			indx = find(iCF==freqs2);
			if ~isempty(indx)
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
			else
				time_last = 0;
			end
		end


		% Shift so all unwrap properly
		start_t = min(time_example);
		t_avg = time_example-start_t;
		% timing = timing - start_t;


		% figure
		% scatter(timing, freqs2, 'filled')
		% hold on
		% xline(0)

		% Average over the F0 period
		% period = 1/F0;
		% t_avg = zeros(length(harms), 1);
		% for iCF = 1:length(harms)
		% 	indx = find(iCF==freqs2);
		% 	time = timing(indx);
		% 	time2 = mod(time, period);
		% 	t_avg(iCF) = mean(time2, 'omitnan');
		% end

		% Plot on previous
		max_f = max(freqs2);
		%plot(ax1, t_avg+0.021, harms(1:max_f)/1000)

		% Order by CF
		% nexttile
		% hold on
		% scatter(t_avg*1000,harms(1:max_f)/1000, 'filled')
		% plot(t_avg*1000,harms(1:max_f)/1000);
		% grid on
		% xlabel('Period, 1/F0 (ms)')
		% set(gca, 'FontSize', 16)
		% title('Averaged Peak Times')
		% ylim(ylimits)

		% Calculate df/dt for each nearest neighbor
		dt = diff(t_avg*1000);
		df = diff(harms(1:max_f)/1000);
		dt_zero_ind = find(dt==0);
		v = df./dt;
		v(dt_zero_ind) = 0;

		% nexttile
		% scatter(v,(harms(1:max_f-1)+F0/2)/1000, 'filled')
		% hold on
		% plot(v,(harms(1:max_f-1)+F0/2)/1000)
		% hold on
		% xlabel('Velocity (kHz/ms)')
		% title('Velocity')
		% set(gca, 'FontSize', 16)
		% xline(0, 'k')
		% max_vel = max(abs(v));
		% xlim([-1*max_vel-0.1 max_vel+0.1])
		% grid on
		% v_harms = harms(1:max_f-1)+F0/2;
		% ylim(ylimits)

		% figure
		% t = tiledlayout(1,1);
		% ax1 = axes(t);
		% ax2 = axes(t);
		% ax1.XAxisLocation = 'bottom';
		% ax2.XAxisLocation = 'top';
		% scatter(ax1,v,(harms(1:end-1)+F0/2)/1000, 'filled')
		% xline(ax1, 0)
		% hold(ax1, "on")
		% plot(ax1, v,(harms(1:end-1)+F0/2)/1000)
		% hold on
		% ax2.Color = 'none';
		% scatter(ax2, t_avg*1000,harms/1000, 'filled')
		% plot(ax2, t_avg*1000,harms/1000);

		%% Save velocities

		%savefile = sprintf('%s_F0_%d_Velocity2.mat', target, round(target_F0));
		%save(fullfile(savepath, savefile), 'v', 'v_harms')
	end
end
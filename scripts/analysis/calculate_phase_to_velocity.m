%% example_chirp.m
clear 

%% Load in exaple 

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

ind = 18; %7;
target_F0 = F0s(ind);
target_file = fullfile(fpath,'waveforms', files{ind});



%% Plot 

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

		fig = figure('Renderer', 'painters', 'Position', [100 100 1200 800]);
		ax = axes('Position',[0.05 0.56 0.28 0.40]);
		pwelch(data_part,[],[],[],fs);
		xlim([0 5])
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
		set(gca,'FontSize',8)

		ax = axes('Position',[0.38 0.56 0.28 0.40]);
		stem(freq,A);
		xlabel('Frequency (Hz)')
		ylabel('Magnitude')
		title('Approx. Mag. Spectrum')
		xlim([0 max_freq])
		set(gca,'FontSize',8)


		ax = axes('Position',[0.38 0.07 0.28 0.40]);
		stem(freq,uwphi);
		xlabel('Frequency (Hz)')
		ylabel('Phase (radians)')
		title('Approx. Phase Spectrum')
		xlim([0 max_freq])
		set(gca,'FontSize',8)

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
	window = round(1/decomp_info(parti).F0_actual*60*1000); %round(decomp_info(parti).F0_actual*6); %600;
	ov = window-10; %round(window*0.9833); %590; % or possibly try minus 10 instead?
	[sg,Ftmp,Ttmp] = spectrogram(pin,window,ov,[],fs,'yaxis');
	Ttmp_part(parti) = {Ttmp};
	Ftmp_part(parti) = {Ftmp};
	spect(parti) = {20*log10(abs(sg))};
	fi = 1;
	while Ftmp(fi) < max_freq
		fi = fi+1;
	end

	% Wait on plotting till the very end
	if parti == 1 % spectrogram of original vowel for comparison
		[sg2,Ftmp2,Ttmp2] = spectrogram(data,window,[],[],fs,'yaxis');
		spect2 = 20*log10(abs(sg2));
		ax = axes('Position',[0.71 0.07 0.28 0.40]);
		spec_image2 = pcolor(1000*Ttmp2,Ftmp2(1:fi),spect2(1:fi,:));
		shading interp
		ylabel('Freq (Hz)')
		xlabel('Time (ms)')
		title('Original Spectrogram')
		set(gca,'FontSize',8)
		xlim([2 1/target_F0*10*1000])

	end

end

for parti = 1:num_parts
	this_spect = spect{parti};
	Ttmp = Ttmp_part{parti};
	Ftmp = Ftmp_part{parti};
	ax = axes('Position',[0.71+(0.13825*(parti-1)) 0.56 0.275/num_parts 0.40]);
	spec_image = pcolor(1000*Ttmp,Ftmp(1:fi),this_spect(1:fi,:));
	shading interp
	if parti == 1
		ylabel('Freq (Hz)')
		xlabel('Time (ms)')
		title(['Synthesized Spectrogram, w=' num2str(window)])
		set(gca,'FontSize',8)
		xlim([2 1/target_F0*10*1000])

	else
		set(gca,'YTickLabel',[]);
	end
end

%% Plot 

% Plot magnitude of harmonics
fig = figure('Renderer', 'painters', 'Position', [70,292,700,500]);
ax = axes('Position',[0.12 0.14 0.30 0.80]);
%subplot(1, 5, 1:2)
stem(freq,20*log10(A),'BaseValue',-inf,'Color','black');
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title(['Mag. Spectrum, F0=' num2str(round(target_F0)) 'Hz'])
xlim([0 max_freq])
ylim([-75 -15])
view(90,90)
set(gca, 'XDir','reverse')
set(gca,'FontSize',16)
set(gca,'YTickLabel',[]);
set(gca,'box','off')

hold on
num_harms = 25;
for iharm = 1:num_harms
	xline(target_F0*iharm, ':');
end

% Plot 'spectrogram' of phase
ax = axes('Position',[0.43 0.14 0.40 0.80]);
%ax = subplot(1, 5, 3:4);
pcolor(1000*Ttmp,Ftmp(1:fi),this_spect(1:fi,:));
shading interp
xlabel('Time (ms)')
set(gca,'YTickLabel',[]);
title('Synthesized Spectrogram')

%alternative caxis settings
% colormap(flipud(gray))
% caxis([-40 50])

% Show 10 periods
xlim([10 1/target_F0*10*1000])
colormap(gray)
ax.CLim = [30, max(max(this_spect(1:fi,:)))];
set(gca,'FontSize',16)
set(gca,'box','off')


% Plot MASD 
f = Ftmp(1:fi);
t = Ttmp;
spec = this_spect(1:fi,:);
spec_diff = diff(spec, 1);
spec_abs = abs(spec_diff);
spec_int = trapz(t, spec_abs, 2);

ax = axes('Position',[0.83 0.14 0.12 0.80]);
plot(f(1:end-1), spec_int, 'k')
set(gca, 'XDir','reverse')
view(90,90)
hold on
num_harms = 25;
for iharm = 1:num_harms
	xline(target_F0*iharm, ':');
end
xticks([])
title('MASD')
set(gca,'FontSize',16)
set(gca,'box','off')
ylim([0 2])


%% Finding peaks using a threshold

ind = 18; %7;
F0 = F0s(ind);
period_lim = (1/F0)*6;
ylimits = [0 4000];
xlimits = [0.02 period_lim+0.02];

% Cut down to avoid onsets/offsets
Fs = fs;
this_spect = this_spect(1:fi,:);
Ftmp = Ftmp(1:fi);
%t_lims = [0.025 0.225]*Fs;
%an = an_sout(:,t_lims(1):t_lims(2)-1);
%t = linspace(0, 0.2, size(an,2));

% Full neurogram
figure('Position',[1840,653,1047,459])
tiledlayout(1, 9)
ax = nexttile([1,3]);
%s = surf(t, CFs, an);
s = surf(Ttmp,Ftmp,this_spect);
s.EdgeColor = 'none';
view(0,90)
xlim(xlimits)
title('AN Neurogram')
xlabel('Time (s)')
ylim(ylimits)
set(gca, 'FontSize', 16)
colorbar
ax.CLim = [30, max(max(this_spect(1:fi,:)))];
colormap(gray)


% Using findpeaks function
nexttile([1,3])
hold on
peaks_all = NaN(size(this_spect,1), size(this_spect,2));
for k1 = 1:size(this_spect,1)
	[pks,loc] = findpeaks(this_spect(k1,:));
	P{k1} = [pks; loc];
	%CF = CFs(k1);

	peaks = NaN(1,size(this_spect, 2));
	peaks(loc) = Ftmp(k1); %CF;
	scatter(t, peaks, 15, 'filled' , 'MarkerEdgeColor','k', 'MarkerFaceColor',...
		'k')
	peaks_all(k1, :) = peaks;
end
xlim(xlimits)
ylim(ylimits)
set(gca, 'FontSize', 16)
title('Peaks using findpeaks')
xlabel('Time (s)')
%yticks([100 200 500 1000 2000 5000 10000])
num_harms = round(4000/F0);
for ii = 1:num_harms 
	yline(F0*ii, 'color', [0.4 0.4 0.4])
end

% Using findpeaks function
nexttile([1,3])
hold on

num_harms = floor(4000/F0);
harms = F0:F0:4000;
peaks_all = NaN(size(this_spect,1), size(this_spect,2));
npts = floor(1/F0*Fs);
for j1 = 1:num_harms

	[~, k1] = min(abs(Ftmp-harms(j1)));
	[pks,loc] = findpeaks(this_spect(k1,:),'MinPeakDistance',25);
	P{k1} = [pks; loc];
	%CF = CFs(k1);

	peaks = NaN(1,size(this_spect, 2));
	peaks(loc) = Ftmp(k1); %CF;
	scatter(t, peaks, 15, 'filled' , 'MarkerEdgeColor','k', 'MarkerFaceColor',...
		'k')
	peaks_all(k1, :) = peaks;
end
xlim(xlimits)
ylim(ylimits)
set(gca, 'FontSize', 16)
title('Peaks at each harmonic')
xlabel('Time (s)')
for ii = 1:num_harms 
	yline(F0*ii, 'color', [0.4 0.4 0.4])
end

%% Mean absolute spatial derivative 

%To extract these cues, we use a spatial derivative operation that simulates
% a hypothetical lateral inhib- itory mechanism operating on the spatiotemporal
% pattern of AN activity (Shamma, 1985). Specifically, we compute the
% point-by-point difference between adjacent rows in Figure 1, and then 
% integrate the absolute value of the difference over time. The resulting 
% “mean absolute spatial deriva- tive” (MASD) shows local maxima at CFs 
% corresponding to the frequen- cies of harmonics 2– 6, whereas the average
% firing rate is mostly saturated at this stimulus level [50 dB sound 
% pressure level (SPL) per component]. Thus, the model predicts that 
% spatiotemporal pitch cues persist at stim- ulus levels where the 
% rate–place representation is degraded.

f = Ftmp(1:fi);
t = 1000*Ttmp;
spec = this_spect(1:fi,:);

% Difference between adjacent rows 
spec_diff = diff(spec, 1);
figure
tiledlayout(1,3)
nexttile
pcolor(1000*Ttmp,f(1:end-1),spec_diff);
shading interp
ylabel('Frequency (Hz)')
xlabel('Time (ms)')
title('Row differences')
colorbar

% Get absolute value of the differece
spec_abs = abs(spec_diff);
nexttile
pcolor(1000*Ttmp,f(1:end-1),spec_abs);
shading interp
ylabel('Frequency (Hz)')
xlabel('Time (ms)')
title('|Row differences|')
colorbar

% Integrate over time 
spec_int = trapz(spec_abs, 2);
nexttile
plot(f(1:end-1), spec_int)


%% Find velocity from the phase plots 
% 
% freq2 = freq/1000;
% 
% figure
% ax = axes('Position',[0.38 0.07 0.28 0.40]);
% stem(freq2,uwphi);
% xlabel('Frequency (Hz)')
% ylabel('Phase (radians)')
% title('Approx. Phase Spectrum')
% xlim([0 max_freq])
% set(gca,'FontSize',8)
% 
% 
% % Convert frequency to angular frequency
% omega = 2 * pi * freq2;
% 
% % Calculate group delay
% group_delay = -gradient(uwphi, omega);
% 
% % Calculate phase velocity (assuming L = 1 for simplicity)
% L = 1;
% phase_velocity = L ./ (uwphi ./ omega);
% 
% % Plot phase velocity vs frequency
% % figure;
% % plot(freq2, phase_velocity);
% % xlabel('Frequency (Hz)');
% % ylabel('Phase Velocity (m/s)');
% % title('Phase Velocity vs Frequency');
% 
% % Calculate and plot group velocity
% group_velocity = L ./ group_delay;
% 
% 
% % Plot 'spectrogram' of phase
% figure;
% tiledlayout(1, 5, 'TileSpacing','tight')
% nexttile([1, 4])
% pcolor(1000*Ttmp,Ftmp(1:fi),this_spect(1:fi,:));
% hold on
% shading interp
% ylabel('Frequency (Hz)')
% xlabel('Time (ms)')
% title('Synthesized Spectrogram')
% 
% %alternative caxis settings
% % colormap(flipud(gray))
% % caxis([-40 50])
% 
% % Show 10 periods
% xlim([10 1/target_F0*10*1000])
% colormap(gray)
% clim([30, max(max(this_spect(1:fi,:)))])
% set(gca,'FontSize',16)
% set(gca,'box','off')
% hold on
% 
% nexttile
% plot(freq2, group_velocity);
% view(90, 90)
% %xlabel('Frequency (Hz)');
% xticklabels([])
% ylabel('Group Velocity (???)');
% title('Instantaneous Velocity');
% set(gca, 'XDir','reverse')
% set(gca,'FontSize',16)
% grid on

%% 
% 
% f = Ftmp(1:fi);
% t = 1000*Ttmp;
% spec = this_spect(1:fi,:);
% 
% figure
% mean_spec = mean(spec, 2);
% plot(f, mean_spec)
% 
% figure
% max_spec = max(spec,[], 2);
% plot(f, max_spec)
% 
% 
% %%
% vel = diff(spec);
% 
% figure
% pcolor(t,f(1:end-1),vel);
% %shading interp
% xlim([10 1/target_F0*10*1000])
% colormap(gray)
% set(gca,'FontSize',16)
% set(gca,'box','off')
% hold on
% 
% figure
% mean_vel = mean(vel, 2);
% plot(f(1:end-1), mean_vel)
% 
% figure
% max_vel = max(vel,[], 2);
% plot(f(1:end-1), max_vel)
% 
% figure
% sum_vel = sum(vel, 2);
% plot(f(1:end-1), sum_vel)
% 
% %% 
% % 
% % 
% % figure
% % tiledlayout(2, 1)
% % nexttile 
% % stem(freq,20*log10(A),'BaseValue',-inf,'Color','black');
% % xlabel('Frequency (Hz)')
% % ylabel('Magnitude (dB)')
% % title(['Mag. Spectrum, F0=' num2str(round(target_F0)) 'Hz'])
% % xlim([0 max_freq])
% % ylim([-75 -15])
% % view(90,90)
% % set(gca,'FontSize',16)
% % set(gca,'YTickLabel',[]);
% % set(gca,'box','off')
% % 
% % 
% % nexttile
% % pcolor(1000*Ttmp,Ftmp(1:fi),this_spect(1:fi,:));
% % shading interp
% % ylabel('Frequency (Hz)')
% % xlabel('Time (ms)')
% % set(gca,'YTickLabel',[]);
% % title('Synthesized Spectrogram')
% % xlim([10 1/target_F0*10*1000])
% % colormap(gray)
% % ax.CLim = [30, max(max(this_spect(1:fi,:)))];
% 
% %%
% 
% %nexttile
% f = Ftmp(1:fi);
% t = Ttmp;
% spec = this_spect(1:fi,:);
% spec(spec<0) = 0;
% 
% figure
% pcolor(t, f, spec)
% shading interp
% colorbar
% 
% figure
% instfreq(spec, f, t)
% 
% 
% % Just looking at mean response
% % mean_t = mean(spec, 2);
% % figure
% % plot(f, mean_t)
% 
% % Incorrect 
% %[ifq, t_if] = instfreq(spec, f, t);
% %vel = diff(ifq);
% 
% 
% %[ifq, t_if] = instfreq(spec', t);
% 
% % figure;
% % [ifq, t_if] = instfreq(spec, f, t);
% % plot(t_if, ifq);
% % view(2);
% % xlabel('Time (s)');
% % ylabel('Frequency (Hz)');
% % title('Instantaneous Frequency');
% % colorbar;

%% Find velocity frequency by frequency
% 
% f = Ftmp(1:fi);
% t2 = 1000*Ttmp;
% spec = this_spect(1:fi,:);
% 
% f_array = target_F0/2:target_F0:4000;
% for isteps = 4:length(f_array)
% 
% 	freq_lo = f_array(isteps-1);
% 	freq_hi = f_array(isteps);
% 	f_ind = find(f<freq_hi & f>freq_lo);
% 
% 	s_part = spec(f_ind,:);
% 	period = 1/target_F0;
% 	period_array = 0:period:dur;
% 
% 	figure
% 	imagesc(t2, f(f_ind), s_part)
% 	hold on
% 	for ii = 1:length(period_array)
% 		xline(period_array(ii), 'LineWidth',2)
% 	end
% 	set(gca, 'YDir', 'normal');
% 	colorbar
% 
% 
% 	%t_ind = 1:length(t2); %find(t<period_array(3) & t>period_array(2));
% 	s_new = s_part; %(:,t_ind);
% 	t_new = t2; %(t_ind);
% 	f_new = f(f_ind);
% 	s_new(s_new < 0) = 0;
% 
% 	% figure
% 	% imagesc(t_new, f_new, s_new)
% 	% hold on
% 	% set(gca, 'YDir', 'normal');
% 
% 	%s_new(s_new<0) = 0;
% 	[ifq, t_if] = instfreq(s_new, f_new, t_new);
% 	t_if = t_if*1000;
% 	ifq = ifq/1000;
% 	figure
% 	plot(t_if, ifq)
% 
% 	dt = diff(t_if);
% 	difq = diff(ifq); %, 1, 2); % Difference along the time axis
% 	velocity = difq ./ dt; % Velocity in Hz/s
% 	vel(isteps) = max(velocity);
% 
% 	% figure
% 	% plot(velocity)
% 
% end
% 
% figure
% plot(f_array(4:end), vel)

%% calculate_velocity_neurogram
clear 

%% Load in natural timbre file names 

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

%% Plot neurograms 

if ismac
	savepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/AN_neurogram';
else
	savepath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\AN_neurogram';
end

for ii = 25

	% Get single stimulus
	target = extractBefore(files{ii}, '.');
	target_F0 = extractBetween(files{ii}, 'ff.', '.wav');
	savefile = sprintf('%s_F0_%s_AN.mat', target, target_F0{1});
	load(fullfile(savepath, savefile),'AN', 'model_params', 'params')

	% Load in stimulus 
	stim = params.stimulus(ii,:);

	% Analyze AN_sout
	an_sout = squeeze(AN.an_sout);
	CFs = AN.CFs;
	t = linspace(0, params.dur, size(an_sout,2));
	F0 = F0s(ii);
	NHN = CFs/F0; % Calculate neural harmonic number (CF/F0) 
	ylimits = [min(CFs)/F0 25]; % max(CFs)/F0];
	period_lim = 1/F0*5;

	% Analyze spectrogram
	Fs = params.Fs;
    y2 = fft(stim);
    m = abs(y2);
    mdB = 20*log10(m);
    f = (0:length(y2)-1)*Fs/length(y2);
    mdB(mdB<0) = 0;
    f(f>max(CFs)) = [];
    mdB = mdB(1:length(f));
	NHN_stim = f/F0;

	% Plot MASD, exclude onset
	t_lims = [0.025 0.225]*Fs;
	an_sout2 = an_sout(:,t_lims(1):t_lims(2)-1);
	t2 = linspace(0, 0.2, size(an_sout2,2));
	spec_diff = diff(an_sout2, 1);
	%spec_abs = abs(spec_diff);
	MASD = trapz(t2, spec_diff, 2);

	% Plot stimulus spectrogram
	figure('Position',[515,275,1082,573])
	tiledlayout(1, 7, 'TileSpacing','none')
	ax1 = nexttile([1,3]);
	plot(NHN_stim, mdB, 'k')
	view(90,90)
	set(gca, 'XDir','reverse')
	title('Stimulus Spectrogram')
	ylabel('Neural Harmonic Number (CF/F0)')
	xlim(ylimits)
	ylabel('Magnitude')
	set(gca, 'FontSize', 16)

	% Plot neurogram 
	ax2 = nexttile([1,3]);
	imagesc(NHN, t, an_sout2')
	set(gca, 'YDir', 'normal');
	set(gca, 'XDir', 'reverse');
	view(90,90)
	ylim([0 period_lim])
	title('AN Neurogram')
	ylabel('Time (s)')
	xlim(ylimits)
	set(gca, 'FontSize', 16)

	% Plot MASD
	ax3 = nexttile;
	plot(NHN(1:end-1), MASD, 'k')
	set(gca, 'XDir','reverse')
	view(90,90)
	hold on
	xticks([])
	title('MASD')
	set(gca,'FontSize',16)
	xlim(ylimits)

	linkaxes([ax1 ax2 ax3],'x')
end

%% Finding peaks using a threshold 
ylimits = [min(CFs) 25*F0]; 
xlimits = [0 period_lim];

% Cut down to avoid onsets/offsets
t_lims = [0.025 0.225]*Fs;
an = an_sout(:,t_lims(1):t_lims(2)-1);
t = linspace(0, 0.2, size(an,2));

% Full neurogram 
figure('Position',[1840,653,1047,459])
tiledlayout(1, 7)
nexttile([1,3])
s = surf(t, CFs, an);
s.EdgeColor = 'none';
view(0,90)
xlim(xlimits)
title('AN Neurogram')
xlabel('Time (s)')
ylim(ylimits)
set(gca, 'FontSize', 16)
colorbar
set(gca, 'YScale', 'log')
yticks([100 200 500 1000 2000 5000 10000])

% Using findpeaks function
nexttile([1,3])
hold on
peaks_all = NaN(size(an,1), size(an,2));
for k1 = 1:size(an,1)
    [pks,loc] = findpeaks(an(k1,:));
    P{k1} = [pks; loc];
	CF = CFs(k1);

	peaks = NaN(1,size(an, 2));
	peaks(loc) = CF;
	scatter(t, peaks, 5, 'filled' , 'MarkerEdgeColor','k', 'MarkerFaceColor',...
		'k')
	peaks_all(k1, :) = peaks;
end
xlim(xlimits)
ylim(ylimits)
set(gca, 'FontSize', 16)
title('Peaks using findpeaks')
set(gca, 'YScale', 'log')
xlabel('Time (s)')
yticks([100 200 500 1000 2000 5000 10000])

nexttile()
peaks_MASD = peaks_all;
peaks_MASD(isnan(peaks_MASD)) = 0;
spec_diff = diff(peaks_MASD, 1);
%spec_abs = abs(spec_diff);
MASD = trapz(t2, spec_diff, 2);
plot(CFs(1:end-1), MASD, 'k')
set(gca, 'XDir','reverse')
view(90,90)
hold on
title('MASD')
set(gca,'FontSize',16)


% %% Plot on top of each other 
% limits = [min(CFs) 25*F0]; 
% 
% % Neurogram 
% figure
% s = surf(t, CFs, an_sout);
% s.EdgeColor = 'none';
% view(0,90)
% xlim([0.05 period_lim+0.05])
% ylim(limits)
% 
% % Using findpeaks function
% hold on
% for k1 = 1:size(an_sout,1)
%     [pks,loc] = findpeaks(an_sout(k1,:));
%     P{k1} = [pks; loc];
% 	CF = CFs(k1);
% 
% 	peaks = NaN(1,36000);
% 	peaks(loc) = CF;
% 	scatter3(t, peaks, 1000, 20, 'filled' , 'MarkerEdgeColor','k', 'MarkerFaceColor',...
% 		'k')
% 	peaks_all(k1, :) = peaks;
% end
% set(gca, 'FontSize', 16)
% title('Peaks using findpeaks')
% set(gca, 'YScale', 'log')

%% Find CF bounds to split neurogram into sections

dur = 0.2;
num_periods = floor(dur/(1/F0));
num_in_row = sum(peaks_all'>0);

% Find boundaries 
CF_diff = diff(num_in_row);
CF_bound_ind = find(CF_diff>num_periods-1);
CF_bound = CFs(CF_bound_ind);

% Plot num in each row 
figure
tiledlayout(1, 2)
nexttile
plot(CFs, num_in_row)
hold on
for ii = 1:length(CF_bound)
	xline(CF_bound(ii))
end
set(gca, 'XScale', 'log')
title('Number of peaks in each CF row')
xlabel('CFs (Hz)')
xticks([100 200 500 1000 2000 5000 10000])
xlim([125 F0*25])

% Split data up into lowest area 
peaks_sub = peaks_all(1:CF_bound_ind(1), :)./1000;

% Plot 
nexttile
scatter(t, peaks_sub,10, 'filled' , 'MarkerEdgeColor','k',...
	'MarkerFaceColor','k')
set(gca, 'YScale', 'log')
xlabel('Time (ms)')
ylabel('CFs (kHz)')
title('Subset')
set(gca, 'FontSize', 16)
xlim(xlimits)

%% Try with multiple periods 
t = linspace(0, 0.2, size(an,2));
%CFs_sub = CFs(1:CF_bound_ind(1));
CFs_sub = CFs(CF_bound_ind(1)+1:CF_bound_ind(2));

% Split data up into target area 
%peaks_sub = peaks_all(1:CF_bound_ind(1), :)./1000;
peaks_sub = peaks_all(CF_bound_ind(1)+1:CF_bound_ind(2), :)./1000;

% Plot 
figure('Position',[612,750,1180,462])
tiledlayout(1, 3)
nexttile
scatter(t, peaks_sub,10, 'filled' , 'MarkerEdgeColor','k',...
	'MarkerFaceColor','k')
set(gca, 'YScale', 'log')
xlabel('Time (ms)')
ylabel('CFs (kHz)')
title('Subset')
set(gca, 'FontSize', 16)
xlim(xlimits)

% Get an array of timing and frequency information
peaks_sub(isnan(peaks_sub)) = 0;
[rows, cols, values] = find(peaks_sub);
timing = t(cols);
freqs = values;
freqs2 = rows;

% Shift so the highest CF is the 'start' of the cycle
[~, start_ind] = max(freqs);
start_t = timing(start_ind) - 0.0001;
timing = timing - start_t;

% Average over the F0 period 
nexttile
period = 1/F0;
for iCF = 1:length(CFs_sub)
	indx = find(iCF==freqs2);
	time = timing(indx);
	time2 = mod(time, period);
	t_avg(iCF) = mean(time2, 'omitnan');
end

% Order by CF
hold on
scatter(t_avg*1000,CFs_sub/1000, 'filled')
plot(t_avg*1000,CFs_sub/1000);
set(gca, 'YScale', 'log')
grid on
xlabel('Period, 1/F0 (ms)')
ylabel('CF (kHz)')
set(gca, 'FontSize', 16)
title('Averaged Peak Times')

% Calculate df/dt for each nearest neighbor 
dt = diff(t_avg*1000);
df = diff(log10(CFs_sub/1000));
dt_zero_ind = find(dt==0);
v = df./dt;
v(dt_zero_ind) = 0;

nexttile
scatter(v,CFs_sub(1:end-1)/1000, 'filled')
hold on
plot(v,CFs_sub(1:end-1)/1000)
hold on
xlabel('Velocity (kHz/ms)')
ylabel('Frequency (kHz)')
title('Velocity')
set(gca, 'FontSize', 16)
xline(0, 'k')
xlim([-0.5 0.5])
grid on
set(gca, 'YScale', 'log')

%ylim([CFs_sub(1) CFs_sub(end)])


%% Try with only one period 
% 
% t_lims = [0.0545 1/F0+0.054];
% [~, t_ind(1)] = min(abs(t-t_lims(1)));
% [~, t_ind(2)] = min(abs(t-t_lims(2)));
% 
% peaks_cut = peaks_all(:,t_ind(1):t_ind(2));
% peaks_sub = peaks_cut(1:CF_bound_ind(1), :);
% CFs_sub = CFs(1:CF_bound_ind)/1000;
% t_sub = t(t_ind(1):t_ind(2))*1000;
% 
% % Get an array of timing and frequency information
% peaks_sub(isnan(peaks_sub)) = 0;
% [~, cols, values] = find(peaks_sub);
% timing = t_sub(cols);
% freqs = values/1000;
% 
% timing(9:10) = [];
% freqs(9:10) = [];
% 
% [freqs, sortIdx] = sort(freqs);
% timing = timing(sortIdx);
% 
% figure
% tiledlayout(1, 2)
% nexttile
% %scatter(timing, freqs, 'filled')
% plot(timing, freqs)
% xlabel('Time (ms)')
% ylabel('CFs (kHz)')
% title('One Period')
% set(gca, 'FontSize', 16)
% grid on
% ylim([CFs_sub(1) CFs_sub(end)])
% 
% % Calculate df/dt for each nearest neighbor 
% dt = diff(timing);
% df = diff(freqs);
% dt_zero_ind = dt==0;
% v = df./dt';
% v(dt_zero_ind) = 0;
% 
% nexttile
% %scatter(v,CFs_sub(1:end-1), 'filled')
% plot(v,CFs_sub(1:end-1))
% hold on
% xlabel('Velocity (kHz/ms)')
% ylabel('Frequency (kHz)')
% title('Velocity')
% set(gca, 'FontSize', 16)
% xline(0, 'k')
% xlim([-0.5 0.5])
% grid on
% ylim([CFs_sub(1) CFs_sub(end)])


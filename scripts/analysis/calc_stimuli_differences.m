%% calc_stimuli_difference
clear

%% Get list of all timbre stimuli (bassoon)

base = getPathsNT;
targets = {'Oboe', 'Bassoon'};
tuning = readtable(fullfile(base, 'Tuning.xlsx')); % Load tuning once
d = dir(fullfile(base, 'waveforms', '*.wav'));
all_files = sort({d.name});
results = struct();
for iinstr = 1:2
	target = targets{iinstr};

	% Find files containing the target instrument name
	is_target = contains(all_files, target);
	target_files = all_files(is_target);

	% Extract note names
	note_names = extractBetween(target_files, 'ff.', '.');

	% Match note names to tuning table and get frequencies
	[~, idx_in_tuning] = ismember(note_names, tuning.Note);
	pitch_order = tuning.Frequency(idx_in_tuning);

	% Sort by frequency
	[~, order] = sort(pitch_order);
	sorted_files = target_files(order);
	sorted_note_names = note_names(order);

	% Store in results struct
	results(iinstr).files = sorted_files;
	results(iinstr).note_names = sorted_note_names;
end

%% Extract overlapping F0s/stimuli numbers (Bb3 to D5)

% 16 overlapping stimuli
n_overlap = 16;
lower_lim = 'Bb3';
high_lim = 'D5';

% Bassoon index
ind_l = find(strcmp(lower_lim, results(2).note_names));
ind_h = find(strcmp(high_lim, results(2).note_names));
bassoon_files = results(2).files(ind_l:ind_h);
F0s = pitch_order(order(ind_l:ind_h));

% Oboe index
ind_l = find(strcmp(lower_lim, results(1).note_names));
ind_h = find(strcmp(high_lim, results(1).note_names));
oboe_files = results(1).files(ind_l:ind_h);
oboe_files(2) = [];


%% Load and plot all / extract envelopes

figure('Position',[68,278,1126,686])
tiledlayout(4, 4, "TileSpacing","compact", 'Padding','compact');
colors = {"#0072BD", "#D95319"};
for j1 = 1:n_overlap

	% Load in data
	[oboe, ~] = audioread(fullfile(base, 'waveforms', oboe_files{j1}));
	[bassoon, Fs] = audioread(fullfile(base, 'waveforms', bassoon_files{j1}));
	name =  extractBetween(oboe_files{j1}, 'ff.','.');

	% Analyze spectrograms & plot
	nexttile
	hold on
	for ii = 1:2

		% Calculate spectra
		F0 = F0s(j1);
		dist = round(F0/4);
		if ii == 1
			stim = oboe;
		else
			stim = bassoon;
		end
		y2 = fft(stim);
		m = abs(y2);
		mdB = 20*log10(m);
		f = (0:length(y2)-1)*Fs/length(y2);
		mdB(mdB<0) = 0;
		f(f>Fs/2) = [];
		mdB = mdB(1:length(f))';

		% Plot
		plot(f/1000, mdB, 'LineWidth', 1.5, 'Color',...
			colors{ii});
		[pks, locs] = findpeaks(mdB, 'MinPeakDistance', dist);
		freqs = f(locs);
		if freqs(1)<F0-10
			plot(freqs(2:end)./1000, pks(2:end), '--', 'LineWidth', 1.5, 'Color',...
				colors{ii});
		else
			plot(freqs./1000, pks, '--', 'LineWidth', 1.5, 'Color',...
				colors{ii});
		end

		envelope(j1).F0 = F0;
		if ii == 1
			envelope(j1).freqs_o = freqs(2:end);
			envelope(j1).peaks_o = pks(2:end);
		else
			envelope(j1).freqs_b = freqs(2:end);
			envelope(j1).peaks_b = pks(2:end);
		end

	end

	% Figure Properties
	xlim([150 8000]./1000)
	xticks([0.2 0.5 1 2 5 10])
	set(gca, 'XScale', 'log')
	hold on
	box on
	ax = gca;
	ax.XGrid='on';
	ax.XMinorGrid='off';

	if ismember(j1, 13:16)
		xlabel('Frequency (kHz)')
	else
		xticklabels([])
	end

	if ismember(j1, [1, 5, 9, 13])
		ylabel('Mag. (dB SPL)')
	else
		yticklabels([])
	end

	if j1 == 1
		legend('Oboe', '','Bassoon')
	end

	text(0.05, 0.95, name{1}, 'Units', 'normalized', ...
		'VerticalAlignment', 'top', 'FontSize',16)
	set(gca,'fontsize',14)
end

%% Calculate differences between envelope

figure('Position',[68,278,1126,686])
hold on
%tiledlayout(4, 4, "TileSpacing","compact", 'Padding','compact');

freq = logspace(log10(50), log10(10000), 100);
env_diff = NaN(16, length(freq));
for ii = 1:16

	% Interpolate envelopes to match
	F0 = envelope(ii).F0;
	
	peaks_o = interp1(envelope(ii).freqs_o, envelope(ii).peaks_o, freq);
	peaks_b = interp1(envelope(ii).freqs_b, envelope(ii).peaks_b, freq);

	% figure
	% plot(freq, peaks_o)
	% hold on
	% plot(freq, peaks_b)

	% Subtract
	env_diff(ii, :) = peaks_o - peaks_b;

	% Plot output
	%nexttile
	plot(freq, env_diff(ii, :))
	hold on
	yline(0, 'k')
	set(gca, 'xscale', 'log')
	ylabel('Oboe Env - Bassoon Env (dB SPL)')
	xlabel('Frequency (Hz)')
end

%% Average envelope difference 

avg_env_diff = mean(env_diff, 'omitnan');

figure
plot(freq, avg_env_diff)
yline(0, 'k')
set(gca, 'xscale', 'log')
grid on

%% Save data 

save(fullfile(base, 'stimuli_avg_env_diff.mat'), 'avg_env_diff', "freq", "env_diff")
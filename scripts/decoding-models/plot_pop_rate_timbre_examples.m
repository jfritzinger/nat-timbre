%%
clear

%% Load in data

filepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/model_comparisons';
load(fullfile(filepath, 'Data_NT.mat'), 'nat_data')

%% Create matrix of stimulus spectra 

if ismac
	fpath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data';
else
	fpath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\';
end
tuning = readtable(fullfile(fpath, 'Tuning.xlsx')); % Load in tuning

% Get bassoon stimulus
target = 'Bassoon';
listing = dir(fullfile(fpath, 'waveforms', ['*' target '*.wav']));
files = {listing.name};
note_names = extractBetween(files, 'ff.', '.');
[~, index] = ismember(note_names, tuning.Note);
F0s = tuning.Frequency(index);
[F0s, order] = sort(F0s);
files = files(order); % Sort files and frequencies by pitch
for ii = 1:size(files, 2)
	[y(ii,:), Fs] = audioread(fullfile(fpath,'waveforms', files{ii}));
end
t = linspace(0, length(y)/Fs, length(y));
for ii = 1:size(files, 2)

	% Calculate spectra
	stim = y(ii,:);
	F0 = F0s(ii);
	dist = round(F0/4);
	y2 = fft(stim);
	m = abs(y2);
	mdB = 20*log10(m);
	f = (0:length(y2)-1)*Fs/length(y2);
	mdB(mdB<0) = 0;
	f(f>Fs/2) = [];
	mdB = mdB(1:length(f))';
	log_f = log10(f(:,2:end));
	mdB_all(ii,:) = mdB;
	f_all(ii, :) = f;


end

% Get oboe stimulus 
target = 'Oboe';
listing = dir(fullfile(fpath, 'waveforms', ['*' target '*.wav']));
files = {listing.name};
note_names = extractBetween(files, 'ff.', '.');
[~, index] = ismember(note_names, tuning.Note);
F0s = tuning.Frequency(index);
[F0s, order] = sort(F0s);
files = files(order); % Sort files and frequencies by pitch
for ii = 1:size(files, 2)
	[y(ii,:), Fs] = audioread(fullfile(fpath,'waveforms', files{ii}));
end
t = linspace(0, length(y)/Fs, length(y));
for ii = 1:size(files, 2)

	% Calculate spectra
	stim = y(ii,:);
	F0 = F0s(ii);
	dist = round(F0/4);
	y2 = fft(stim);
	m = abs(y2);
	mdB = 20*log10(m);
	f = (0:length(y2)-1)*Fs/length(y2);
	mdB(mdB<0) = 0;
	f(f>Fs/2) = [];
	mdB = mdB(1:length(f))';
	log_f = log10(f(:,2:end));
	mdB_all_oboe(ii,:) = mdB;
	f_all_oboe(ii, :) = f;

end

%% Sort data 

sesh = [];
for ii = 1:length(nat_data)
	rate = nat_data(ii).bass_rate;
	if ~isempty(rate)
		sesh = [sesh ii];
	end
end
num_data = length(sesh);
bass_rate = [nat_data(sesh).bass_rate];
bass_CFs = [nat_data(sesh).CF];

sesh = [];
for ii = 1:length(nat_data)
	rate = nat_data(ii).oboe_rate;
	if ~isempty(rate)
		sesh = [sesh ii];
	end
end
num_data = length(sesh);
oboe_rate = [nat_data(sesh).oboe_rate];
oboe_CFs = [nat_data(sesh).CF];

%% Plot population of vector strengths 

ind_b = 25:40;
ind_o = [1 3:17];
for ind = 1:16

	figure
	tiledlayout(3, 1)
	F0 = F0s(ind);
	dist = round(F0/4);

	nexttile
	hold on
	plot(f_all(ind_b(ind), :)/1000, mdB_all(ind_b(ind),:), 'LineWidth', 1.5, 'Color',"#0072BD");
	[pks, locs] = findpeaks(mdB_all(ind_b(ind),:), 'MinPeakDistance', dist);
	freqs = f(locs);
	plot(freqs/1000, pks, '--', 'LineWidth', 1.5, 'Color',"#0072BD");

	plot(f_all_oboe(ind_o(ind), :)/1000, mdB_all_oboe(ind_o(ind),:), 'LineWidth', 1.5, 'Color',"#D95319");
	[pks, locs] = findpeaks(mdB_all_oboe(ind_o(ind),:), 'MinPeakDistance', dist);
	freqs = f(locs);
	plot(freqs/1000, pks, '--', 'LineWidth', 1.5, 'Color',"#D95319");

	ylabel('Mag. (dB SPL)')
	xlim([0.1 10])
	xlabel('Frequency (kHz)')
	set(gca, 'xscale', 'log')
	xticks([0.1 0.2 0.5 1 2 5 10])

	VS = bass_rate;
	CFs = bass_CFs;
	VS2 = VS(ind_b(ind),:);
	nexttile 
	scatter(CFs/1000, VS2, 'filled', 'MarkerEdgeColor','k')
	set(gca, 'xscale', 'log')
	xticks([0.1 0.2 0.5 1 2 5 10])
	xlim([0.1 10])
	ylabel('Norm Rate')
	%ylim([0 1])
	xlabel('CF (kHz)')
	grid on

	VS = oboe_rate;
	CFs = oboe_CFs;
	VS2 = VS(ind_o(ind),:);
	nexttile
	scatter(CFs/1000, VS2, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceColor',	"#D95319")
	set(gca, 'xscale', 'log')
	xticks([0.1 0.2 0.5 1 2 5 10])
	xlim([0.1 10])
	ylabel('Norm Rate')
	%ylim([0 1])
	xlabel('CF (kHz)')
	grid on
end


%% timbre_pop_example
clear 
save_fig = 1;

%% Load in data 

[base, datapath, savepath, ppi] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')
load(fullfile(base, 'model_comparisons', 'Pop_Rate_Timbre_All.mat'), ...
	"pop_rate_timbre")
load(fullfile(base, 'model_comparisons','Neuron_Rate_Timbre_All.mat'),...
	"neuron_rate_timbre")
accuracy_rate = [neuron_rate_timbre.accuracy]*100;

sesh = [];
for ii = 1:length(nat_data)
	rate = nat_data(ii).bass_rate;
	rate2 = nat_data(ii).oboe_rate;
	if ~isempty(rate) && ~isempty(rate2)
		sesh = [sesh ii];
	end
end
num_data = length(sesh);

beta_weights = pop_rate_timbre.trainedClassifier.ClassificationSVM.Beta;
%beta_weights = pop_rate_timbre.imp.ImportanceMean;

%% Create figure 

figure('Position',[50, 50, 4.2*ppi, 4.2*ppi])
tiledlayout(3, 1, "TileSpacing","none", 'Padding','compact')
linewidth = 2;
fontsize = 18;
titlesize = 20;
legsize = 20;
scattersize = 40;
colorsData = {"#0072BD", "#D95319"};
colorsMTF = {'#648FFF', '#DC267F', '#785EF0', '#FFB000'};
colorsTimbre = '#1b9e77';


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
for ind = 8

	F0 = F0s(ind);
	dist = round(F0/4);

	nexttile
	hold on
	plot(f_all(ind_b(ind), :)/1000, mdB_all(ind_b(ind),:), 'LineWidth', 1, 'Color',"#0072BD");
	[pks, locs] = findpeaks(mdB_all(ind_b(ind),:), 'MinPeakDistance', dist);
	freqs = f(locs);
	plot(freqs/1000, pks, '--', 'LineWidth', 1.5, 'Color',"#0072BD");

	plot(f_all_oboe(ind_o(ind), :)/1000, mdB_all_oboe(ind_o(ind),:), 'LineWidth', 1, 'Color',"#D95319");
	[pks, locs] = findpeaks(mdB_all_oboe(ind_o(ind),:), 'MinPeakDistance', dist);
	freqs = f(locs);
	plot(freqs/1000, pks, '--', 'LineWidth', 1, 'Color',"#D95319");

	ylabel('Mag.')
	xlim([0.1 10])
	xlabel('Frequency (kHz)')
	set(gca, 'xscale', 'log')
	xticks([0.1 0.2 0.5 1 2 5 10])
	xticklabels([])
	set(gca, 'fontsize', fontsize)
	grid on

	VS = bass_rate;
	CFs = bass_CFs;
	VS2 = VS(ind_b(ind),:);
nexttile 
	scatter(CFs/1000, VS2, scattersize, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.5)
	set(gca, 'xscale', 'log')
	xticks([0.1 0.2 0.5 1 2 5 10])
	xlim([0.1 10])
	xticklabels([])
	yticks([0 100])
	grid on
	set(gca, 'fontsize', fontsize)
	hleg = legend('Bassoon', 'fontsize', legsize, 'location', 'northwest', 'box', 'off');
	hleg.ItemTokenSize = [8, 8];

	VS = oboe_rate;
	CFs = oboe_CFs;
	VS2 = VS(ind_o(ind),:);
nexttile
scatter(CFs/1000, VS2, scattersize, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceColor',	"#D95319", 'MarkerFaceAlpha',0.5)
	set(gca, 'xscale', 'log')
	xticks([0.1 0.2 0.5 1 2 5 10])
	xtickangle(0)
	xlim([0.1 10])
	yticks([0 100])
	ylabel('        Rate (sp/s)')
	%ylim([0 1])
	xlabel('CF (kHz)')
	grid on
	set(gca, 'fontsize', fontsize)
	hleg = legend('Oboe', 'fontsize', legsize, 'location', 'northwest', 'box', 'off');
	hleg.ItemTokenSize = [8, 8];
end

%% Annotate and position

%% Save figure 

if save_fig == 1
	filename = 'timbre_pop_examples';
	save_figure_MARC(filename)
end
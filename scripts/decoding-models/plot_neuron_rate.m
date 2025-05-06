%% plot_neuron_rate_F0
clear

%% Set paths

[base, ~, savepath, ppi] = getPathsNT();

%% Load and plot results for F0

figure('Position',[560,563,780,285])
tiledlayout(1, 4)
linewidth = 2;
instruments = {'Bassoon', 'Oboe'};
for ii = 1:2

	% Load in data
	filename = ['Neuron_Rate_F0_' instruments{ii}];
	filepath = fullfile(base, 'model_comparisons',filename);
	load(filepath, "neuron_rate_F0")

	% Plot accuracy of each neuron
	nexttile
	accuracy = [neuron_rate_F0.accuracy]*100;
	CFs = [neuron_rate_F0.CF];
	edges = linspace(0, 100, 101);
	histogram(accuracy,edges)
	hold on

	% Plot chance line 
	chance = 1/length(neuron_rate_F0(1).rate)*100;
	xline(chance, 'k', 'LineWidth',linewidth)

	% Mean and median 
	xline(mean(accuracy), 'r', 'LineWidth',linewidth)
	xline(median(accuracy), 'r--', 'LineWidth',linewidth)

	% Best 
	scatter(max(accuracy), 0, 'filled')

	% Labels 
	legend('', sprintf('Chance = %.2f%%', chance), ...
		sprintf('Mean = %.2f%%', mean(accuracy)), ...
		sprintf('Median = %.2f%%', median(accuracy)),...
		sprintf('Best = %.2f%%', max(accuracy)))
	ylabel('# Neurons')
	xlabel('Prediction Accuracy (%)')
	title([instruments{ii} ' F0 Predictions'])
	xlim([0 25])
	grid on

	% Plot accuracy vs CF 
	nexttile
	scatter(CFs, accuracy, 'filled')
	set(gca, 'xscale', 'log')
	xlabel('CFs (Hz)')
	ylabel('Accuracy')

end

%% Load and plot results for timbre separate 

% Load in data 
filename = 'Neuron_Rate_Timbre_Separate.mat';
filepath = fullfile(base, 'model_comparisons',filename);
load(filepath, "neuron_rate_timbre")

% Load in timbre F0s
pitch = getBassoonF0s();

% Arrange data
accuracy = NaN(180, 16);
for itarget = 1:16
	accuracy(:,itarget) = [neuron_rate_timbre(:,itarget).accuracy]*100;
end
CFs = [neuron_rate_timbre(:,1).CF];

% Plot accuracy of each neuron
figure('Position',[1872,432,1190,580])
tiledlayout(4, 4, 'TileSpacing','tight')
linewidth = 1;
for ii = 1:16

	% Plot accuracy
	nexttile
	histogram(accuracy(:,ii),21)
	hold on

	% Plot chance line 
	xline(50, 'k', 'LineWidth',linewidth)

	% Mean and median 
	xline(mean(accuracy(:,ii)), 'r', 'LineWidth',linewidth)
	xline(median(accuracy(:,ii)), 'r--', 'LineWidth',linewidth)

	% Best 
	scatter(max(accuracy(:,ii)), 0, 'filled')

	hleg = legend('', sprintf('Chance = 50%%'), ...
		sprintf('Mean = %.2f%%', mean(accuracy(:,ii))), ...
		sprintf('Median = %.2f%%', median(accuracy(:,ii))),...
		sprintf('Best = %.2f%%', max(accuracy(:,ii))), ...
		'Location','northwest');
	hleg.ItemTokenSize = [10, 10];
	title(sprintf('F0 = %0.0f Hz', pitch(neuron_rate_timbre(1, ii).ind_b)))
	xlim([0 100])
	if ismember(ii, [1, 5, 9, 13])
		ylabel('# Neurons')
	end
	if ii >= 13
		xlabel('Accuracy (%)')
	end
end

% Plot accuracy of each neuron vs CF
figure('Position',[1872,432,1190,580])
tiledlayout(4, 4, 'TileSpacing','tight')
linewidth = 1;
for ii = 1:16

	% Plot accuracy
	nexttile
	scatter(CFs, accuracy(:,ii),"filled", 'MarkerEdgeColor','k')
	hold on
	title(sprintf('F0 = %0.0f Hz', pitch(neuron_rate_timbre(1, ii).ind_b)))
	if ismember(ii, [1, 5, 9, 13])
		ylabel('Accuracy (%)')
	end
	ylim([0 100])
	set(gca, 'xscale', 'log')
	if ii >= 13
		xlabel('CF (Hz)')
	end
end

%% Load and plot results for timbre all F0s together 

% Load in data 
filename = 'Neuron_Rate_Timbre_All.mat';
filepath = fullfile(base, 'model_comparisons',filename);
load(filepath, "neuron_rate_timbre")

% Figure
figure('Position',[560,585,1023,263])
tiledlayout(1, 2)
linewidth = 2;

nexttile
accuracy = [neuron_rate_timbre.accuracy]*100;
CFs = [neuron_rate_timbre.CF];
edges = linspace(0, 100, 101);
histogram(accuracy,edges)
hold on

% Plot chance line
xline(50, 'k', 'LineWidth',linewidth)

% Mean and median
xline(mean(accuracy), 'r', 'LineWidth',linewidth)
xline(median(accuracy), 'r--', 'LineWidth',linewidth)

% Best
scatter(max(accuracy), 0, 'filled')

% Labels
legend('', sprintf('Chance = 50%%'), ...
	sprintf('Mean = %.2f%%', mean(accuracy)), ...
	sprintf('Median = %.2f%%', median(accuracy)),...
	sprintf('Best = %.2f%%', max(accuracy)))
ylabel('# Neurons')
xlabel('Prediction Accuracy (%)')
title('Instrument Predictions')
xlim([0 100])
grid on

% Plot accuracy vs CF
nexttile
scatter(CFs, accuracy, 'filled')
set(gca, 'xscale', 'log')
xlabel('CFs (Hz)')
ylabel('Accuracy')


%% FUNCTIONS 


function pitch = getBassoonF0s()
target = 'Bassoon';
if ismac
	fpath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data';
else
	fpath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\';
end
tuning = readtable(fullfile(fpath, 'Tuning.xlsx')); % Load in tuning
listing = dir(fullfile(fpath, 'waveforms', '*.wav'));
target_WAV = arrayfun(@(n) contains(listing(n).name, target), 1:numel(listing), 'UniformOutput', false);
wav_nums =  find(cell2mat(target_WAV));
d = dir(fullfile(fpath,'waveforms', '*.wav'));
all_files = sort({d.name});
nfiles = length(wav_nums);

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
pitch = pitch_order(order);

end
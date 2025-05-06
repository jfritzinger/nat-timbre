%% plot_neuron_comparison_timbre
clear 

[base, ~, savepath, ppi] = getPathsNT();


%% 

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
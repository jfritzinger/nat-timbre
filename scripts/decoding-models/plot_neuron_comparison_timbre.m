%% plot_neuron_comparison_timbre
clear 

[base, ~, ~, ~] = getPathsNT();

%% 

% Load in data 
filename = 'Neuron_Rate_Timbre_All.mat';
filepath = fullfile(base, 'model_comparisons',filename);
load(filepath, "neuron_rate_timbre")

% Figure
figure('Position',[560,585,1023,500])
tiledlayout(2, 2)
linewidth = 2;

nexttile
accuracy_rate = [neuron_rate_timbre.accuracy]*100;
CFs = [neuron_rate_timbre.CF];
edges = linspace(0, 100, 101);
histogram(accuracy_rate,edges)
hold on

% Plot chance line
xline(50, 'k', 'LineWidth',linewidth)

% Mean and median
xline(mean(accuracy_rate), 'r', 'LineWidth',linewidth)
xline(median(accuracy_rate), 'r--', 'LineWidth',linewidth)

% Best
scatter(max(accuracy_rate), 0, 'filled')

% Labels
legend('', sprintf('Chance = 50%%'), ...
	sprintf('Mean = %.2f%%', mean(accuracy_rate)), ...
	sprintf('Median = %.2f%%', median(accuracy_rate)),...
	sprintf('Best = %.2f%%', max(accuracy_rate)))
ylabel('# Neurons')
xlabel('Prediction Accuracy (%)')
title('Instrument Predictions')
xlim([0 100])
grid on

% Plot accuracy vs CF
nexttile
scatter(CFs, accuracy_rate, 'filled')
set(gca, 'xscale', 'log')
xlabel('CFs (Hz)')
ylabel('Accuracy')

%% 

% Load in data 
[base, datapath, savepath, ppi] = getPathsNT();
filename = 'Neuron_Time_Timbre_All.mat';
filepath = fullfile(base, 'model_comparisons',filename);
load(filepath, "neuron_time_timbre")

% Figure
nexttile
accuracy_time = [neuron_time_timbre.accuracy]*100;
CFs = [neuron_time_timbre.CF];
edges = linspace(0, 100, 101);
histogram(accuracy_time,edges)
hold on

% Plot chance line
xline(50, 'k', 'LineWidth',linewidth)

% Mean and median
xline(mean(accuracy_time), 'r', 'LineWidth',linewidth)
xline(median(accuracy_time), 'r--', 'LineWidth',linewidth)

% Best
scatter(max(accuracy_time), 0, 'filled')

% Labels
legend('', sprintf('Chance = 50%%'), ...
	sprintf('Mean = %.2f%%', mean(accuracy_time)), ...
	sprintf('Median = %.2f%%', median(accuracy_time)),...
	sprintf('Best = %.2f%%', max(accuracy_time)))
ylabel('# Neurons')
xlabel('Prediction Accuracy (%)')
title('Instrument Predictions')
xlim([0 100])
grid on

% Plot accuracy vs CF
nexttile
scatter(CFs, accuracy_time, 'filled')
set(gca, 'xscale', 'log')
xlabel('CFs (Hz)')
ylabel('Accuracy')

%% Look at single neuron differences in predictions 

figure
hold on
scatter(accuracy_rate, accuracy_time, 'filled', 'MarkerEdgeColor','k', ...
	MarkerFaceAlpha=0.5)
plot([0, 100], [0, 100], 'k')
xlim([40 85])
ylim([40 85])
axis square

mdl = fitlm(accuracy_rate, accuracy_time);
x = linspace(0, 100, 20);
y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
plot(x, y, 'r')
legend('Neuron', 'Unity', 'y = 0.48*x+30.47, p=0.0000')
xlabel('Accuracy of Rate Prediction (%)')
ylabel('Accuracy of Timing Prediction (%)')
title('Rate vs Timing Comparison')

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
%% plot_neuron_time
clear

%% Load in data
target = 'Oboe';

pitch = getF0s(target);
msg{1} = sprintf('%0.0f - %0.0f Hz F0', pitch(1), pitch(20));
msg{2} = sprintf('%0.0f - %0.0f Hz F0', pitch(21), pitch(end));

[base, datapath, savepath, ppi] = getPathsNT();
load(fullfile(base, 'model_comparisons', ['Neuron_Time_F0_' target '.mat']), ...
	"neuron_time_F0")

%% Plot accuracy of each neuron
figure('Position',[560,618,798,230])
tiledlayout(1, 3)
linewidth = 1;

accuracy(1,:) = [neuron_time_F0.accuracy]*100;
accuracy(2,:) = [neuron_time_F0.accuracy_low]*100;
accuracy(3,:) = [neuron_time_F0.accuracy_high]*100;

% Plot overall histogram
nexttile
edges = linspace(0, 100, 30);
histogram(accuracy(1,:),edges)
hold on
chance = 1/size(neuron_time_F0(1).C, 1)*100;
xline(chance, 'k', 'LineWidth',linewidth)
xline(mean(accuracy(1,:)), 'r', 'LineWidth',2)
scatter(max(accuracy(1,:)), 0, 'filled')
ylabel('# Neurons')
xlabel('Prediction Accuracy (%)')
title([target ' Prediction of F0'])
xlim([0 100])
legend('', sprintf('Chance = %.2f%%', chance), ...
	sprintf('Mean = %.2f%%', mean(accuracy(1,:))), ...
	sprintf('Best = %.2f%%', max(accuracy(1,:))))

% Split histogram into two sections
nexttile
edges = linspace(0, 100, 30);
histogram(accuracy(2,:),edges, 'FaceAlpha',0.5, 'FaceColor','b')
hold on
histogram(accuracy(3,:),edges, 'FaceAlpha',0.5, 'FaceColor','r')
xline(chance, 'k', 'LineWidth',linewidth)
xline(mean(accuracy(2,:)), 'b', 'LineWidth',1)
xline(mean(accuracy(3,:)), 'r', 'LineWidth',1)
ylabel('# Neurons')
xlabel('Prediction Accuracy (%)')
title('Low vs high F0 Predictions')
xlim([0 100])
legend(msg{1}, msg{2}, ...
	sprintf('Chance = %.2f%%', chance), ...
	sprintf('Mean = %.2f%%', mean(accuracy(2,:))), ...
	sprintf('Mean = %.2f%%', mean(accuracy(3,:))))


% Find indices of the best neurons 
[temp,originalpos] = sort(accuracy(1,:), 'descend' );
n = temp(1:3);
best_ind=originalpos(1:3);

% Plot confusion matrix of best neuron
putatives = neuron_time_F0(best_ind(1)).putative;
nexttile
confusionchart(neuron_time_F0(best_ind(1)).C)
title(sprintf('Best prediction, %s, %0.02f%%', putatives, temp(1)))

%% Test accuracy vs neuron characteristics, like CF and MTF 




%% Load and plot results for timbre all F0s together 

% Load in data 
[base, datapath, savepath, ppi] = getPathsNT();
filename = 'Neuron_Time_Timbre_All.mat';
filepath = fullfile(base, 'model_comparisons',filename);
load(filepath, "neuron_time_timbre")

% Figure
figure('Position',[560,585,1023,263])
tiledlayout(1, 2)
linewidth = 2;

nexttile
accuracy = [neuron_time_timbre.accuracy]*100;
CFs = [neuron_time_timbre.CF];
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

function pitch = getF0s(target)

[base, ~, ~, ~] = getPathsNT();
tuning = readtable(fullfile(base, 'Tuning.xlsx')); % Load in tuning
listing = dir(fullfile(base, 'waveforms', '*.wav'));
target_WAV = arrayfun(@(n) contains(listing(n).name, target), 1:numel(listing), 'UniformOutput', false);
wav_nums =  find(cell2mat(target_WAV));
d = dir(fullfile(base,'waveforms', '*.wav'));
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
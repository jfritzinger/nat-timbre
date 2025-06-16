%% pitch_neuron_F0_image
clear
save_fig = 1;

%% Load in data
target = 'Bassoon';

[base, ~, ~, ppi] = getPathsNT();
load(fullfile(base, 'model_comparisons', ['Neuron_Time_F0_' target '.mat']), ...
	"neuron_time_F0")

% load(fullfile(base, 'model_comparisons', 'Neuron_Time_Timbre_All.mat'),...
% 	"neuron_time_timbre")
% neuron_time_F0 = neuron_time_timbre;

%% Plot accuracy of each neuron


figure('Position',[50,50,3.6*ppi,9.4*ppi])
tiledlayout(1, 3)
linewidth = 1;
linewidth = 2;
fontsize = 18;
titlesize = 20;
legsize = 20;
scattersize = 40;


% Plot imagesc for diag for each neuron 

accuracy= [neuron_time_F0.accuracy]*100;
[~,best_ind] = sort(abs(accuracy), 'ascend' );
for ii = 1:287
	C_acc(ii,:) = diag(neuron_time_F0(best_ind(ii)).C);
end

x = getF0s(target);
y = 1:287;

pcolor(x, y, C_acc, 'EdgeColor','none', 'EdgeAlpha',0)
set(gca, 'xscale', 'log')
xticks([60 100 200 350 550])
xlabel('F0 (Hz)')
ylabel('All Neurons')
c = colorbar;
title('F0 predictions\newlinefor all neurons', 'HorizontalAlignment','center')
set(gca, 'fontsize', fontsize)
c.Label.String = '# Accurate Predictions';

% Annotations 



% Save figure 

if save_fig == 1
	filename = 'pich_neuron_F0_image';
	save_figure_MARC(filename)
end

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
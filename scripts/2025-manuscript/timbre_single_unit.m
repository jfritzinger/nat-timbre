%% timbre_single_unit
% This script plots figure 2 of the natural timbre manuscript, which is an
% analysis of decoding models fit to individual neurons
clear
save_fig = 0;

%% Load in data

[base, ~, ~, ppi] = getPathsNT();

% Load in rate models
filepath = fullfile(base, 'model_comparisons','Neuron_Rate_Timbre_All.mat');
load(filepath, "neuron_rate_timbre")

% Load in timing models
filepath_timing = fullfile(base, 'model_comparisons', 'Neuron_Time_Timbre_All.mat');
load(filepath_timing, "neuron_time_timbre")

% Load in all data
load(fullfile(fullfile(base, 'model_comparisons', 'Data_NT_3.mat')), 'nat_data')

%% Set up figure

figure('Position',[50, 50, 6.7*ppi, 3*ppi])
linewidth = 1;
fontsize = 8;
labelsize = 12;
legsize = 6;

%% A, B. Plot examples of accurate and inaccurate


pitch = getBassoonF0s();
for igood = 1:2

	accuracy_rate = [neuron_rate_timbre.accuracy]*100;
	if igood == 1
		[ac, originalpos] = sort(accuracy_rate, 'descend');
		h_ind = [1, 7];
	else
		[ac, originalpos] = sort(accuracy_rate, 'ascend');
		h_ind = [13 19];
	end
	ind_high=originalpos(1:2);
	putative = {neuron_rate_timbre(ind_high).putative};
	ind_b = 25:40;
	ind_o = [1 3:17];
	for iexample = 1:2
		if igood == 1
			h(iexample) = subplot(4, 6, h_ind(iexample));
		else
			h(2+iexample) = subplot(4, 6, h_ind(iexample));
		end

		% Get rates
		pitches = pitch(ind_b);
		index = find(strcmp(putative{iexample}, {nat_data.putative}));
		oboe_rate = nat_data(index).oboe_rate(ind_o);
		oboe_std = nat_data(index).oboe_rate_std(ind_o);
		bass_rate = nat_data(index).bass_rate(ind_b);
		bass_std = nat_data(index).bass_rate_std(ind_b);

		errorbar(pitches, oboe_rate, oboe_std/sqrt(20), 'b', 'LineWidth', linewidth, ...
			'CapSize',3);
		hold on
		errorbar(pitches, bass_rate, bass_std/sqrt(20), 'r', 'LineWidth', linewidth,...
			'CapSize',3);
		box off 
		grid on
		xticks([110 220 440 660])
		yticks([0 25 50 75 100])

		if iexample == 1
			xticklabels([])
			ylabel('Avg Rate (sp/s)               ')
		else
			xlabel('F0 (Hz)')
			xticklabels([110 220 440 660])
		end
		set(gca, 'fontsize', fontsize)

		if iexample == 1 && igood == 2
			hleg = legend('Oboe', 'Bassoon', 'fontsize', legsize);
			hleg.ItemTokenSize = [8, 8];
		end
	end
end


%% C. Plot rate accuracy vs CF (add MTF type)

CFs = [neuron_rate_timbre.CF];
MTFs = {neuron_rate_timbre.MTF};
isMTF(1,:) = strcmp(MTFs, 'BE');
isMTF(2,:) = strcmp(MTFs, 'BS');
isMTF(3,:) = contains(MTFs, 'H');
isMTF(4,:) = strcmp(MTFs, 'F');

h(5) = subplot(4, 6, [2, 3, 8, 9]);
hold on
for iMTF = 1:4
	scatter(CFs(isMTF(iMTF, :)), accuracy_rate(isMTF(iMTF, :)), 10, 'filled', ...
		'MarkerEdgeColor','k')
end
ylim([40 100])
set(gca, 'xscale', 'log')
xticks([100 200 500 1000 2000 5000 10000])
xticklabels([])
ylabel('Rate Accuracy (%)')
set(gca, 'fontsize', fontsize)
grid on

%% D. Plot timing accuracy vs CF (add MTF type)

accuracy_time = [neuron_time_timbre.accuracy]*100;
CFs = [neuron_time_timbre.CF];

h(6) = subplot(4, 6, [14, 15, 20, 21]);
hold on
for iMTF = 1:4
	scatter(CFs(isMTF(iMTF, :)), accuracy_time(isMTF(iMTF, :)), 10, 'filled', ...
		'MarkerEdgeColor','k')
end
hleg = legend('BE', 'BS', 'Hybrid', 'Flat', 'fontsize', legsize);
hleg.ItemTokenSize = [8, 8];
set(gca, 'xscale', 'log')
xlabel('CFs (Hz)')
ylabel('Timing Accuracy (%)')
set(gca, 'fontsize', fontsize)
xticks([100 200 500 1000 2000 5000 10000])
xticklabels([100 200 500 1000 2000 5000 10000]/1000)
grid on
ylim([40 100])

%% E. Plot histogram and scatters of rate vs timing

% Histogram rate and histogram timing
for iplot = 1:2
	if iplot == 1
		h(7) = subplot(4, 6, [4 10 16]);
		accuracy = accuracy_rate;
	else
		h(8) = subplot(4, 6, [23, 24]);
		accuracy = accuracy_time;
	end
	edges = linspace(0, 100, 51);
	if iplot == 1
		histogram(accuracy,edges, 'Orientation','horizontal')
		ylim([40 85])
		ylabel('Timing Prediction Accuracy (%)')
		yline(50, 'k', 'LineWidth',linewidth)
		yline(mean(accuracy), 'r', 'LineWidth',linewidth)
		yline(median(accuracy), 'r--', 'LineWidth',linewidth)
	else
		histogram(accuracy,edges)
		xlim([40 85])
		xlabel('Rate Prediction Accuracy (%)')
		xline(50, 'k', 'LineWidth',linewidth)
		xline(mean(accuracy), 'r', 'LineWidth',linewidth)
		xline(median(accuracy), 'r--', 'LineWidth',linewidth)
	end
	hold on

	%legend('', sprintf('Chance = 50%%'), ...
	%	sprintf('Mean = %.2f%%', mean(accuracy)), ...
	%	sprintf('Median = %.2f%%', median(accuracy)))
	grid on
	set(gca, 'fontsize', fontsize)

end

% Scatter
h(9) = subplot(4, 6, [5, 6, 11, 12, 17, 18]);
hold on
scatter(accuracy_rate, accuracy_time,10, 'filled', 'MarkerEdgeColor','k', ...
	MarkerFaceAlpha=0.5)
plot([0, 100], [0, 100], 'k')
xlim([40 85])
ylim([40 85])

mdl = fitlm(accuracy_rate, accuracy_time);
x = linspace(0, 100, 20);
y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
plot(x, y, 'r')
yticks(0:5:100)
xticks(0:5:100)
hleg = legend('Neuron', 'Unity', ...
	sprintf('y = %0.2f*x+%0.2f, p=%0.04f', mdl.Coefficients{2, 1}, ...
	mdl.Coefficients{1, 1},mdl.Coefficients{2,4}), 'fontsize', legsize);
hleg.ItemTokenSize = [8, 8];
title('Rate vs Timing Comparison')
set(gca, 'fontsize', fontsize)
xticklabels([])
yticklabels([])
grid on


%% Arrange figure and add annotations

left = [0.07 0.3];
bottom = [0.12 0.29 0.59 0.77];
width = 0.15;
height = 0.15;

set(h(1), 'Position', [left(1) bottom(4) width height])
set(h(2), 'Position', [left(1) bottom(3) width height])
set(h(3), 'Position', [left(1) bottom(2) width height])
set(h(4), 'Position', [left(1) bottom(1) width height])

set(h(5), "Position", [left(2), 0.56, 0.23 0.36])
set(h(6), "Position", [left(2), bottom(1), 0.23 0.36])

fig_position = [0.68,0.26,0.29,0.63];
nb_position = [fig_position(1),fig_position(2)-0.14,fig_position(3),0.11];
wb_position = [fig_position(1)-0.07,fig_position(2),0.06,fig_position(4)];
set(h(9), 'Position', fig_position)
set(h(8), 'Position', nb_position)
set(h(7), 'Position', wb_position)

%% Create textbox

labelleft= [0.0, 0.24, 0.55];
labelbottom = [0.48 0.95];
annotation('textbox',[labelleft(1) labelbottom(2) 0.071 0.058],...
	'String','A','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(1) labelbottom(1) 0.071 0.058],...
	'String','B','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(2) labelbottom(2) 0.071 0.058],...
	'String','C','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(3) labelbottom(2) 0.071 0.058],...
	'String','D','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');

%% Save figure 

if save_fig == 1
	filename = 'fig2_timbre_single_unit';
	save_figure(filename)
end


%%

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
index = zeros(1, nfiles);
note_names = extractBetween(files, 'ff.','.');
for ii = 1:nfiles % Find index of each note in tuning spreadsheet
	index(ii) = find(strcmp(note_names(ii), tuning.Note));
end
pitch_order = tuning.Frequency(index); % Get freqs of each note
[~, order] = sort(pitch_order); % Sort freqs
pitch = pitch_order(order);

end
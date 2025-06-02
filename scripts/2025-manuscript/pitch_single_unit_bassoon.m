%% plot_pitch_singleunit
clear
save_fig = 1;

%% Load in data 

target = 'Bassoon';
pitch = getF0s(target);
[base, ~, ~, ppi] = getPathsNT();
load(fullfile(base, 'model_comparisons', ['Neuron_Time_F0_' target '.mat']), ...
	"neuron_time_F0")
load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')

accuracy = [neuron_time_F0.accuracy]*100;
% accuracy(2,:) = [neuron_time_F0.accuracy_low]*100;
% accuracy(3,:) = [neuron_time_F0.accuracy_high]*100;

%% Set up figure 

figure('position', [50 50 6.7*ppi 3*ppi])
scattersize = 5;
titlesize = 9;
fontsize = 8;
labelsize = 12;
legsize = 7;
linewidth = 1;

%% A. 

% Get best single-unit timbre neurons 
filepath = fullfile(base, 'model_comparisons','Neuron_Rate_F0_Bassoon.mat');
load(filepath, "neuron_rate_F0")
accuracy_rate = [neuron_rate_F0.accuracy]*100;
accuracy_time = accuracy;

h(1) = subplot(3, 3, 1);
set(gca, 'fontsize', fontsize)
edges = linspace(0, 20, 21);
histogram(accuracy_rate,edges)
xlabel('Rate Accuracy')
grid on
box off 
set(gca, 'fontsize', fontsize)

h(2) = subplot(3, 3, 2);
set(gca, 'fontsize', fontsize)
edges = linspace(0, 60, 21);
histogram(accuracy_time, edges, 'Orientation','horizontal')
ylabel('Timing Accuracy')
grid on
box off
set(gca, 'fontsize', fontsize)

h(3) = subplot(3, 3, 3);
hold on
scatter(accuracy_rate, accuracy_time, scattersize, 'filled', 'MarkerEdgeColor','k')
plot([0 100], [0 100], 'k')
grid on
xlim([0 20])
ylim([0 60])
set(gca, 'fontsize', fontsize)
xticklabels([])
yticklabels([])

% mdl = fitlm(accuracy_rate, accuracy_time);
% x = linspace(0, 60, 10);
% y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
% plot(x, y, 'r')
% hleg = legend('Neuron', 'p=7.1671e-25');
% hleg.ItemTokenSize = [8, 8];


%% B. CF vs accuracy
h(4) = subplot(3, 3, 4);

CFs = [neuron_time_F0.CF];
scatter(CFs, accuracy(1,:), scattersize, 'filled', 'MarkerEdgeColor','k');
set(gca, 'xscale', 'log')
hold on
set(gca, 'fontsize', fontsize)
xticks([200 500 1000 2000 5000 10000])
xticklabels([0.2 0.5 1 2 5 10])
xlabel('CF (kHz)')
ylabel('Accuracy')

%% C. MTF vs accuracy 
h(5) = subplot(3, 3, 5);

MTFs = {neuron_time_F0.MTF};
MTF_types = unique(MTFs);
hold on
for iMTF = 1:5
	ind = strcmp(MTFs, MTF_types{iMTF});
	accur = accuracy(1,ind);
	num_units = length(accur);

	swarmchart(ones(num_units, 1)*iMTF, accur, scattersize)

	mean_vals(iMTF) = mean(accur);
	std_vals(iMTF) = std(accur)/sqrt(length(accur));
end
errorbar(1:5, mean_vals, std_vals, 'k')
xticks(1:5)
xticklabels(MTF_types)
xlabel('MTF Groups')
ylabel('Accuracy')

% tableMTF = table(MTFs', accuracy(1,:)');
% anova(tableMTF, 'Var2')
% [~,~,stats] = anova1(accuracy(1,:), MTFs);
% [c,~,~,gnames] = multcompare(stats);
grid on
set(gca, 'fontsize', fontsize)

%% D. PCA2 vs accuracy 
h(6) = subplot(3, 3, 6);

VS_all = zeros(3,1);
for ii = 1:297
	VS_1 = nat_data(ii).bass_VS;
	if ~isempty(VS_1)
		PC2_score(ii,:) = nat_data(ii).RVF_PC2;
	else
		PC2_score(ii,:) = 0;
	end
end
PC2_score(PC2_score==0) = [];


% Plot PCA2 score vs beta weights 
scatter(PC2_score, accuracy(1,:), scattersize, 'filled', 'MarkerEdgeColor','k')
xlabel('PC2 RVF score')
ylabel('Accuracy')
hold on

mdl = fitlm(PC2_score, accuracy(1,:));
x = linspace(-2, 1, 20);
y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
plot(x, y, 'r')
hleg = legend('Neuron', ...
	sprintf('p=%0.04f',mdl.Coefficients{2,4}), 'fontsize', 7, ...
	'location', 'northwest', 'box', 'off');
hleg.ItemTokenSize = [8, 8];
grid on
set(gca, 'fontsize', fontsize)


%% E. VS vs accuracy
h(7) = subplot(3, 3, 7);

load(fullfile(base,'model_comparisons', 'Data_NT_3.mat'), 'nat_data')

% Get subset of units 
VS_all = zeros(3,1);
for ii = 1:297
	VS_1 = nat_data(ii).bass_VS;
	if ~isempty(VS_1)
		VS_all(ii,:) = mean(VS_1);
	end
end
VS_all(VS_all==0) = [];

scatter(VS_all, accuracy, scattersize, 'filled', 'MarkerEdgeColor','k');
hold on

mdl = fitlm(VS_all, accuracy);
x = linspace(0, 0.6, 20);
y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
plot(x, y, 'r')
hleg = legend('Neuron', ...
	sprintf('p=%0.04f', mdl.Coefficients{2,4}), 'box', 'off');
hleg.ItemTokenSize = [8, 8];
ylabel('Accuracy (%)')
xlabel('Vector Strength')
set(gca, 'fontsize', fontsize)
ylim([0 62])

%% G. Low vs high 
h(8) = subplot(3, 3, 8);

accuraclo = [neuron_time_F0.accuracy_low]*100;
accuracyhi = [neuron_time_F0.accuracy_high]*100;
chance = 2.5;
edges = linspace(0, 100, 30);
histogram(accuraclo,edges, 'FaceAlpha',0.5, 'FaceColor','b')
hold on
histogram(accuracyhi,edges, 'FaceAlpha',0.5, 'FaceColor','r')
xline(chance, 'k', 'LineWidth',linewidth)
xline(mean(accuraclo), 'b', 'LineWidth',1)
xline(mean(accuracyhi), 'r', 'LineWidth',1)
ylabel('# Neurons')
xlabel('Prediction Accuracy (%)')
xlim([0 100])
% hleg = legend(msg{1}, msg{2}, ...
% 	sprintf('Chance = %.2f%%', chance), ...
% 	sprintf('Mean = %.2f%%', mean(accuraclo)), ...
% 	sprintf('Mean = %.2f%%', mean(accuracyhi)), 'fontsize', legsize);
% hleg.ItemTokenSize = [8, 8];
box off 
set(gca, 'fontsize', fontsize)


%% Arrange plots 

left = linspace(0.42, 0.83, 3);
bottom = linspace(0.13, 0.65, 2);
height = 0.3;
width = 0.15;

fig_position = [0.12,0.23,0.23,0.67];
nb_position = [fig_position(1),fig_position(2)-0.14,fig_position(3),0.11];
wb_position = [fig_position(1)-0.07,fig_position(2),0.06,fig_position(4)];
set(h(3), 'Position', fig_position)
set(h(1), 'Position', nb_position)
set(h(2), 'Position', wb_position)

% set(h(3), 'Position', [left(1) bottom(3) width height/3])
% set(h(2), 'Position', [left(1) bottom(3)+height/3 width height/3])
% set(h(1), 'Position', [left(1) bottom(3)+height/3*2 width height/3])

set(h(4), 'position', [left(1) bottom(2) width height])
set(h(5), 'position', [left(2) bottom(2) width height])
set(h(6), 'position', [left(3) bottom(2) width height])
set(h(7), 'position', [left(1) bottom(1) width height])
set(h(8), 'position', [left(2) bottom(1) width height])

%% Set labels 

labelleft= [0.0 linspace(0.365, 0.775, 3)];
labelbottom = [0.46 0.95];
annotation('textbox',[labelleft(1) labelbottom(2) 0.071 0.058],...
	'String','A','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(2) labelbottom(2) 0.071 0.058],...
	'String','B','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(3) labelbottom(2) 0.071 0.058],...
	'String','C','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(4) labelbottom(2) 0.071 0.058],...
	'String','D','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(2) labelbottom(1) 0.071 0.058],...
	'String','E','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(3) labelbottom(1) 0.071 0.058],...
	'String','F','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');

%% Save figure 

if save_fig == 1
	filename = 'fig4_pitch_single_unit_bassoon';
	save_figure(filename)
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
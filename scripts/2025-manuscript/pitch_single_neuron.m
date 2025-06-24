%%pitch_single_unit 
clear
save_fig = 0;

%% Load in data 

[base, ~, ~, ppi] = getPathsNT();

% Load in oboe
load(fullfile(base, 'model_comparisons', 'Neuron_Time_F0_Oboe.mat'), ...
	"neuron_time_F0")
neuron_time_F0_oboe = neuron_time_F0;
load(fullfile(base, 'model_comparisons','Neuron_Rate_F0_Oboe.mat'), ...
	"neuron_rate_F0")
neuron_rate_F0_oboe = neuron_rate_F0;

% Load in bassoon
load(fullfile(base, 'model_comparisons', 'Neuron_Time_F0_Bassoon.mat'), ...
	"neuron_time_F0")
load(fullfile(base, 'model_comparisons','Neuron_Rate_F0_Bassoon.mat'), ...
	"neuron_rate_F0")

accuracy_rate = [neuron_rate_F0.accuracy]*100;
accuracy_time = [neuron_time_F0.accuracy]*100;
accuracy_rate_oboe = [neuron_rate_F0_oboe.accuracy]*100;
accuracy_time_oboe = [neuron_time_F0_oboe.accuracy]*100;
CFs = [neuron_time_F0.CF];
MTFs = {neuron_time_F0.MTF};

%% Create figure

figure('Position',[560,563,6.9*ppi,6*ppi])
linewidth = 1;
fontsize = 8;
legsize = 7;
labelsize = 12; 
colorsData = {"#0072BD", "#D95319"};
colorsMTF = {'#648FFF', '#DC267F', '#785EF0', '#FFB000'};
colorsPitch = '#0072B2';
scattersize = 10;

%% A. Bassoon rate vs timing 

h(1) = subplot(4, 3, 1);
set(gca, 'fontsize', fontsize)
edges = linspace(0, 20, 21);
histogram(accuracy_rate,edges, 'FaceColor', colorsPitch)
xlabel('Rate Accuracy')
grid on
box off 
set(gca, 'fontsize', fontsize)

h(2) = subplot(4, 3, 2);
set(gca, 'fontsize', fontsize)
edges = linspace(0, 60, 21);
histogram(accuracy_time, edges, 'Orientation','horizontal', 'FaceColor', colorsPitch)
ylabel('Timing Accuracy')
grid on
box off
set(gca, 'fontsize', fontsize)

h(3) = subplot(4, 3, 3);
hold on
scatter(accuracy_rate, accuracy_time, scattersize, 'filled',...
	'MarkerEdgeColor','k', 'MarkerFaceColor', colorsPitch, 'MarkerFaceAlpha',0.5)
plot([0 100], [0 100], 'k')
grid on
xlim([0 20])
ylim([0 60])
set(gca, 'fontsize', fontsize)
xticklabels([])
yticklabels([])
title('Bassoon')



%% B. Oboe rate vs timing 

for itype = 1:2

	h(3+itype) = subplot(4, 3, 3+itype);
	edges = linspace(0, 100, 101);
	chance = 1/length(neuron_rate_F0_oboe(1).rate)*100;
	hold on
	if itype == 2
		accuracy = accuracy_time_oboe;
		histogram(accuracy,edges, 'Orientation','horizontal')
		yline(chance, 'k', 'LineWidth',linewidth)
		yline(mean(accuracy), 'r', 'LineWidth',linewidth)
		ylim([0 25])
		ylabel('Timing Prediction Accuracy (%)')
	else
		accuracy = accuracy_rate_oboe;
		histogram(accuracy,edges)
		xline(chance, 'k', 'LineWidth',linewidth)
		xline(mean(accuracy), 'r', 'LineWidth',linewidth)
		xlim([0 25])
		xlabel('Rate Prediction Accuracy (%)')
	end
	grid on
	set(gca, 'fontsize', fontsize)
end


h(6) = subplot(4, 3, 6);
scatter(accuracy_rate, accuracy_time, 10, 'filled', 'MarkerEdgeColor','k', ...
	'MarkerFaceAlpha',0.5)
hold on
plot([0, 100], [0, 100], 'k')
xlim([0 25])
ylim([0 25])

% mdl = fitlm(accuracy_rate, accuracy_time);
% x = linspace(0, 100, 20);
% y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
% plot(x, y, 'r')
yticks(0:5:100)
xticks(0:5:100)
% hleg = legend('Neuron', 'Unity', ...
% 	sprintf('y = %0.2f*x+%0.2f, p=%0.04f', mdl.Coefficients{2, 1}, ...
% 	mdl.Coefficients{1, 1},mdl.Coefficients{2,4}), 'fontsize', legsize);
% hleg.ItemTokenSize = [8, 8];
title('Oboe')
set(gca, 'fontsize', fontsize)
xticklabels([])
yticklabels([])
grid on


%% C. CFs

h(7) = subplot(4, 3, 7);

scatter(CFs, accuracy_time, scattersize, 'filled', 'MarkerEdgeColor',...
	'k', 'MarkerFaceColor', colorsPitch, 'MarkerFaceAlpha',0.5);
set(gca, 'xscale', 'log')
hold on
set(gca, 'fontsize', fontsize)
xticks([200 500 1000 2000 5000 10000])
xticklabels([0.2 0.5 1 2 5 10])
xlabel('CF (kHz)')
ylabel('Bassoon Accuracy (%)')
grid on

%% D. MTFs

h(8) = subplot(4, 3, 8);

MTF_types = unique(MTFs);
hold on
for iMTF = 1:4
	if iMTF == 4
		ind = strcmp(MTFs, MTF_types{iMTF}) | strcmp(MTFs, MTF_types{iMTF+1});
	else
		ind = strcmp(MTFs, MTF_types{iMTF});
	end
	accur = accuracy_time(ind);
	num_units = length(accur);

	RGB = hex2rgb(colorsMTF{iMTF});
	swarmchart(ones(num_units, 1)*iMTF, accur, scattersize, RGB)

	mean_vals(iMTF) = mean(accur);
	std_vals(iMTF) = std(accur)/sqrt(length(accur));
end
errorbar(1:4, mean_vals, std_vals, 'k')
xticks(1:4)
xticklabels({'BE', 'BS', 'F', 'H'})
xlabel('MTF Groups')
ylabel('Bassoon Accuracy (%)')

% tableMTF = table(MTFs', accuracy(1,:)');
% anova(tableMTF, 'Var2')
% [~,~,stats] = anova1(accuracy(1,:), MTFs);
% [c,~,~,gnames] = multcompare(stats);
grid on
set(gca, 'fontsize', fontsize)


%% E. Vector strength

h(9) = subplot(4, 3, 9);
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

scatter(VS_all, accuracy_time, scattersize, 'filled',...
	'MarkerEdgeColor','k', 'MarkerFaceColor', colorsPitch, 'MarkerFaceAlpha',0.5);
hold on

mdl = fitlm(VS_all, accuracy_time);
x = linspace(0, 0.6, 20);
y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
plot(x, y, ':k')
% hleg = legend('Neuron', ...
% 	sprintf('p=%0.04f', mdl.Coefficients{2,4}), 'box', 'off');
% hleg.ItemTokenSize = [8, 8];
ylabel('Bassoon Accuracy (%)')
xlabel('Vector Strength')
set(gca, 'fontsize', fontsize)
ylim([0 62])
grid on

%% F. Imagesc 

h(10) = subplot(4, 3, 10);

[~,best_ind] = sort(abs(accuracy_time), 'ascend' );
for ii = 1:287
	C_acc(ii,:) = diag(neuron_time_F0(best_ind(ii)).C);
end

x = getF0s('Bassoon');
y = 1:287;

pcolor(x, y, C_acc, 'EdgeColor','none', 'EdgeAlpha',0)
set(gca, 'xscale', 'log')
xticks([60 100 200 350 550])
xlabel('F0 (Hz)')
ylabel('All Neurons')
c = colorbar;
title('Bassoon F0 predictions\newlinefor all neurons',...
	'HorizontalAlignment','center')
set(gca, 'fontsize', fontsize)
c.Label.String = '# Accurate Predictions';

%% Arrange 

left = [0.12 0.48 0.74];
bottom = linspace(0.08, 0.73, 3);
height = 0.22;
width = 0.18;

fig_position = [left(1),0.65,0.28,0.3];
nb_position = [fig_position(1),fig_position(2)-0.08,fig_position(3),0.06];
wb_position = [fig_position(1)-0.07,fig_position(2),0.06,fig_position(4)];
set(h(3), 'Position', fig_position)
set(h(1), 'Position', nb_position)
set(h(2), 'Position', wb_position)

fig_position = [left(1),bottom(1)+0.08,0.28,0.3];
nb_position = [fig_position(1),fig_position(2)-0.08,fig_position(3),0.06];
wb_position = [fig_position(1)-0.07,fig_position(2),0.06,fig_position(4)];
set(h(6), 'Position', fig_position)
set(h(4), 'Position', nb_position)
set(h(5), 'Position', wb_position)


set(h(7), 'position', [left(2) bottom(3) width height])
set(h(8), 'position', [left(2) bottom(2) width height])
set(h(9), 'position', [left(2) bottom(1) width height])
set(h(10), 'position', [left(3) bottom(1) width 0.83])

%% Set labels 

labelleft= [0 0.42 0.7];
labelbottom = linspace(0.29, 0.94, 3);
annotation('textbox',[labelleft(1) labelbottom(3) 0.071 0.058],...
	'String','A','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(1) 0.47 0.071 0.058],...
	'String','B','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(2) labelbottom(3) 0.071 0.058],...
	'String','C','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(2) labelbottom(2) 0.071 0.058],...
	'String','D','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(2) labelbottom(1) 0.071 0.058],...
	'String','E','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(3) labelbottom(3) 0.071 0.058],...
	'String','F','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');

%% Others!

figure
[~,best_ind] = sort(abs(accuracy_time_oboe), 'ascend' );
C_acc = zeros(length(best_ind), 35);
for ii = 1:length(best_ind)
	C_acc(ii,:) = diag(neuron_time_F0_oboe(best_ind(ii)).C);
end

x = getF0s('Oboe');
y = 1:length(best_ind);

pcolor(x, y, C_acc, 'EdgeColor','none', 'EdgeAlpha',0)
set(gca, 'xscale', 'log')
xticks([60 100 200 350 550])
xlabel('F0 (Hz)')
ylabel('All Neurons')
c = colorbar;
title('Oboe F0 predictions TIME\newlinefor all neurons',...
	'HorizontalAlignment','center')
set(gca, 'fontsize', fontsize)
c.Label.String = '# Accurate Predictions';


%% Save figure 

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

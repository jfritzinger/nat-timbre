%% timbre_rate_timing_comparison
clear
save_fig = 1;

%% Load in data

[base, ~, ~, ppi] = getPathsNT();

% Load in rate models
filepath = fullfile(base, 'model_comparisons','Neuron_Rate_Timbre_All.mat');
load(filepath, "neuron_rate_timbre")

% Load in timing models
filepath_timing = fullfile(base, 'model_comparisons', 'Neuron_Time_Timbre_All.mat');
load(filepath_timing, "neuron_time_timbre")

%% Set up figure

accuracy_rate = [neuron_rate_timbre.accuracy];
accuracy_time = [neuron_time_timbre.accuracy];


figure('Position',[50, 50, 5.5*ppi, 5.5*ppi])
linewidth = 2;
fontsize = 18;
titlesize = 20;
legsize = 20;
scattersize = 40;
colorsData = {"#0072BD", "#D95319"};
colorsMTF = {'#648FFF', '#DC267F', '#785EF0', '#FFB000'};
colorsTimbre = '#1b9e77';

% Histogram rate and histogram timing
for iplot = 1:2
	if iplot == 1
		h(7) = subplot(4, 6, [4 10 16]);
		accuracy = accuracy_rate;
	else
		h(8) = subplot(4, 6, [23, 24]);
		accuracy = accuracy_time;
	end
	edges = linspace(0, 1, 51);
	if iplot == 1
		histogram(accuracy,edges, 'Orientation','horizontal', 'FaceColor', colorsTimbre)
		ylim([0.40 0.85])
		ylabel('Timing Accuracy')
		yline(0.50, 'color', [0.4 0.4 0.4], 'LineWidth',linewidth)
		yticks(0:0.05:1)
		xticks([0 0.30])
		%yline(mean(accuracy), 'r', 'LineWidth',linewidth)
		%yline(median(accuracy), 'r--', 'LineWidth',linewidth)
	else
		histogram(accuracy,edges, 'FaceColor', colorsTimbre)
		xlim([0.40 0.85])
		xlabel('Rate Accuracy')
		xline(0.50, 'color', [0.4 0.4 0.4], 'LineWidth',linewidth)
		xticks(0:0.05:1)
		yticks([0 0.30])
		%xline(mean(accuracy), 'r', 'LineWidth',linewidth)
		%xline(median(accuracy), 'r--', 'LineWidth',linewidth)
	end
	box off
	hold on
	grid on
	set(gca, 'fontsize', fontsize)

end

% Scatter
h(9) = subplot(4, 6, [5, 6, 11, 12, 17, 18]);
hold on
scatter(accuracy_rate, accuracy_time,scattersize, 'filled', 'MarkerEdgeColor','k', ...
	'MarkerFaceAlpha', 0.5, 'MarkerFaceColor', colorsTimbre)
plot([0, 1], [0, 1], 'k', 'linewidth', 2)
plot([0 0.50], [0.50, 0.50], 'color', [0.4 0.4 0.4], 'linewidth', 2)
plot([0.50 0.50], [0, 0.50], 'color', [0.4 0.4 0.4], 'linewidth', 2)
xlim([0.40 0.85])
ylim([0.40 0.85])

mdl = fitlm(accuracy_rate, accuracy_time);
x = linspace(0, 1, 20);
y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
plot(x, y, '--k', 'linewidth', 2)
yticks(0:0.05:1)
xticks(0:0.05:1)
hleg = legend('', '','Chance','', ...
	sprintf('p=%0.04f',mdl.Coefficients{2,4}), 'fontsize', legsize, ...
	'position', [0.693666698736137,0.338541370738636,0.283817951959545,0.1487], 'box', 'off');
hleg.ItemTokenSize = [20, 20];
title('Rate vs Timing Comparison')
set(gca, 'fontsize', fontsize)
xticklabels([])
yticklabels([])
grid on


% Arrange figure and add annotations
fig_position = [0.3,0.33,0.68,0.6];
nb_position = [fig_position(1),fig_position(2)-0.14,fig_position(3),0.11];
wb_position = [fig_position(1)-0.14,fig_position(2),0.12,fig_position(4)];
set(h(9), 'Position', fig_position)
set(h(8), 'Position', nb_position)
set(h(7), 'Position', wb_position)

%% Save figure 

if save_fig == 1
	filename = 'timbre_rate_timing_comparison';
	save_figure_MARC(filename)
end
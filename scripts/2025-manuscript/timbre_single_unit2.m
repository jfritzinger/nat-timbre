% timbre_single_unit2

clear
save_fig = 0;

%% Load in data

[base, ~, ~, ppi] = getPathsNT();

% Load in rate models
filepath = fullfile(base, 'model_comparisons','Neuron_Rate_Timbre_All.mat');
load(filepath, "neuron_rate_timbre")
accuracy_rate = [neuron_rate_timbre.accuracy];

% Load in timing models
filepath_timing = fullfile(base, 'model_comparisons', 'Neuron_Time_Timbre_All.mat');
load(filepath_timing, "neuron_time_timbre")
accuracy_time = [neuron_time_timbre.accuracy];

% Load in all data
load(fullfile(fullfile(base, 'model_comparisons', 'Data_NT_3.mat')), 'nat_data')

%% Set up figure

figure('Position',[50, 50, 6.7*ppi, 3*ppi])
linewidth = 1;
fontsize = 8;
labelsize = 12;
legsize = 6;
colorsData = {"#0072BD", "#D95319"};
colorsMTF = {'#648FFF', '#DC267F', '#785EF0', '#FFB000'};
colorsTimbre = '#1b9e77';
scattersize = 12;


%% A. Rate vs timing 

for iplot = 1:2
	if iplot == 1
		h(1) = subplot(2, 3, 1);
		accuracy = accuracy_rate;
	else
		h(2) = subplot(2, 3, 2);
		accuracy = accuracy_time;
	end
	edges = linspace(0, 1, 51);
	if iplot == 1
		histogram(accuracy,edges, 'Orientation','horizontal', 'FaceColor', colorsTimbre)
		ylim([0.40 0.85])
		ylabel('Timing Prediction Accuracy')
		yline(0.50, 'color', [0.4 0.4 0.4], 'LineWidth',linewidth)
		yticks(0:0.05:1)
		xticks([0 30])
		%yline(mean(accuracy), 'r', 'LineWidth',linewidth)
		%yline(median(accuracy), 'r--', 'LineWidth',linewidth)
	else
		histogram(accuracy,edges, 'FaceColor', colorsTimbre)
		xlim([0.40 0.85])
		xlabel('Rate Prediction Accuracy')
		xline(0.50, 'color', [0.4 0.4 0.4], 'LineWidth',linewidth)
		xticks(0:0.05:1)
		yticks([0 30])
		%xline(mean(accuracy), 'r', 'LineWidth',linewidth)
		%xline(median(accuracy), 'r--', 'LineWidth',linewidth)
	end
	box off
	hold on

	%legend('', sprintf('Chance = 50%%'), ...
	%	sprintf('Mean = %.2f%%', mean(accuracy)), ...
	%	sprintf('Median = %.2f%%', median(accuracy)))
	grid on
	set(gca, 'fontsize', fontsize)

end

% Scatter
h(3) = subplot(2, 3, 3);
hold on
scatter(accuracy_rate, accuracy_time, scattersize, 'filled', 'MarkerEdgeColor','k', ...
	'MarkerFaceAlpha', 0.5, 'MarkerFaceColor', colorsTimbre)
plot([0, 1], [0, 1], 'k')
plot([0 0.50], [0.50, 0.50], 'color', [0.4 0.4 0.4])
plot([0.50 0.50], [0, 0.50], 'color', [0.4 0.4 0.4])
xlim([0.40 0.85])
ylim([0.40 0.85])

mdl = fitlm(accuracy_rate, accuracy_time);
x = linspace(0, 1, 20);
y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
plot(x, y, ':k')
yticks(0:0.05:1)
xticks(0:0.05:1)
hleg = legend('Neuron', 'Unity','Chance','', ...
	sprintf('y = %0.2f*x+%0.2f, p=%0.04f', mdl.Coefficients{2, 1}, ...
	mdl.Coefficients{1, 1},mdl.Coefficients{2,4}), 'fontsize', legsize, ...
	'position', [0.26,0.28,0.1877,0.1487], 'box', 'off');
hleg.ItemTokenSize = [8, 8];
title('Rate vs Timing Comparison')
set(gca, 'fontsize', fontsize)
xticklabels([])
yticklabels([])
grid on

%% B. CF distributions

CFs = [neuron_rate_timbre.CF];
MTFs = {neuron_rate_timbre.MTF};
MTF_types = unique(MTFs);

h(4) = subplot(2, 3, 4);
scatter(CFs, accuracy_time, scattersize, 'filled', 'MarkerEdgeColor','k', ...
	'MarkerFaceColor',colorsTimbre, 'MarkerFaceAlpha',0.5);
hold on
set(gca, 'xscale', 'log')

mdl = fitlm(CFs, accuracy_time);
x = linspace(300, 14000, 50);
y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
plot(x, y, ':k')
hleg = legend('Neuron', ...
	sprintf('p=%0.04f',mdl.Coefficients{2,4}), 'fontsize', legsize, ...
	'location', 'northwest', 'box', 'off');
hleg.ItemTokenSize = [8, 8];
ylabel('Accuracy')
xlabel('CFs')
xticks([200 500 1000 2000 5000 10000])
xticklabels([0.2 0.5 1 2 5 10])
set(gca, 'fontsize', fontsize)
grid on
ylim([0.38 1])

%% C. MTF Distributions

h(5) = subplot(2, 3, 5);
[weights_ordered, order_ind] = sort(accuracy_time);
hold on
for iMTF = 1:4
	if iMTF == 4
		ind = strcmp(MTFs, MTF_types{iMTF}) | strcmp(MTFs, MTF_types{iMTF+1});
	else
		ind = strcmp(MTFs, MTF_types{iMTF});
	end
	[weights_ordered, order_ind] = sort(accuracy_time(ind));
	num_units = length(weights_ordered);

	RGB = hex2rgb(colorsMTF{iMTF});
	swarmchart(ones(num_units, 1)*iMTF, weights_ordered, scattersize, RGB, 'filled', 'MarkerFaceAlpha',0.5)

	mean_vals(iMTF) = mean(weights_ordered);
	std_vals(iMTF) = std(weights_ordered)/sqrt(length(weights_ordered));
end
errorbar(1:4, mean_vals, std_vals, 'k')
xticks(1:4)
xticklabels({'BE', 'BS', 'F', 'H'})
xlabel('MTF Groups')
ylabel('Accuracy')
set(gca, 'fontsize', fontsize)
grid on
ylim([0.38 1])

%% D. Bin size comparisons 

load(fullfile(base, 'model_comparisons', 'Neuron_Time_Timbre_All_Shuffled.mat'), ...
	"neuron_time_timbre")
acc_shuffled = [neuron_time_timbre.accuracy];
load(fullfile(base, 'model_comparisons', 'Neuron_Time_Timbre_All_Coarse150.mat'), ...
	"neuron_time_timbre")
acc_150ms = [neuron_time_timbre.accuracy];
acc_025ms = accuracy_time;

h(6) = subplot(2, 3, 6);
accur_all = [acc_shuffled; acc_025ms; acc_150ms; accuracy_rate];
boxplot(accur_all')
hold on
swarmchart(ones(length(acc_shuffled), 1), acc_shuffled, scattersize, 'filled', 'MarkerFaceAlpha',0.4)
swarmchart(ones(length(acc_025ms), 1)*2, acc_025ms, scattersize, 'filled', 'MarkerFaceAlpha',0.4)
swarmchart(ones(length(acc_150ms), 1)*3, acc_150ms, scattersize, 'filled', 'MarkerFaceAlpha',0.4)
swarmchart(ones(length(accuracy_rate), 1)*4, accuracy_rate, scattersize, 'filled', 'MarkerFaceAlpha',0.4)

xlabel('Groups')
ylabel('Accuracy')
xticks(1:4)
xticklabels({'Shuffle', '0.25 ms', '150 ms', 'Rate'})
box off 
set(gca, 'fontsize', fontsize)

% [p12, ~] = ranksum(acc_shuffled, acc_025ms);  
% [p13, ~] = ranksum(acc_shuffled, acc_150ms);  
% [p23, ~] = ranksum(acc_025ms, acc_100ms);  
% adjusted_p = [p12, p13, p23] * 3; % Bonferroni adjustment
% 
% [p, tbl, stats] = anova1(accur_all');
% [c, m, h, gnames] = multcompare(stats, 'CType', 'hsd');


%% Arrange 

left = [0.07 0.52 0.78];
bottom = [0.14 0.6];
width = 0.17;
height = 0.3;

fig_position = [0.16,0.28,0.29,0.63];
nb_position = [fig_position(1),fig_position(2)-0.14,fig_position(3),0.11];
wb_position = [fig_position(1)-0.07,fig_position(2),0.06,fig_position(4)];
set(h(3), 'Position', fig_position)
set(h(2), 'Position', nb_position)
set(h(1), 'Position', wb_position)

set(h(4), 'Position', [left(2) bottom(2) width height])
set(h(5), "Position", [left(2), bottom(1), width height])
set(h(6), "Position", [left(3), bottom(1)+0.05, width 0.7])

%% Annotate 


%% Save 
% fig 4.4 timbre_singleunit_timing
clear 
% save_fig = 1;

%% Load in data 

[base, ~, ~, ppi] = getPathsNT();
load(fullfile(base, 'model_comparisons_revised', 'Neuron_Time_Timbre_Distances_L1.mat'), ...
	"neuron_time_timbre")
load(fullfile(base,'model_comparisons_revised', 'Data_NT_3.mat'), 'nat_data')


num_data = size(neuron_time_timbre, 2);
min_dis = [0.1 0.25 0.5 0.75 1 1.5 2 2.5 3 3.5 4 4.5 5 10 15 20 30 40 50 100 150 300];
accuracy_all = zeros(length(min_dis), num_data);
for i = 1:num_data
	accuracy_all(:,i) = [neuron_time_timbre(:,i).accuracy];
end

%% Set up figure

figure('Position',[50 50 6*ppi, 4.7*ppi])
scattersize = 15;
titlesize = 9;
fontsize = 8;
labelsize = 12;
legsize = 7;
linewidth = 1;
spont_color = [0.4 0.4 0.4];
CF_color = [0.7 0.7 0.7];
colorsTimbre = '#1b9e77';

%% Stats for all 

h(1) = subplot(3, 4, 1);
accuracy_box = accuracy_all([1, 2, 3, 4, 5, 7, 9, 11, 13:22],:);

boxplot(accuracy_box', 'BoxStyle','filled', 'OutlierSize',2, 'MedianStyle','target')
%xticks(1:22)
xticklabels(min_dis([1, 2, 3, 4, 5, 7, 9, 11, 13:22]))
xlabel('PSTH Bin Size (ms)')
ylabel('Accuracy')
grid on
set(gca, 'fontsize', fontsize)
box off 

%% Example PSTHs 

colors = {"#0072BD", "#D95319"};
accuracy_time = accuracy_all(16,:);
[acc_sorted, index] = sort(accuracy_time, 'descend');
iii = [1 2 215 210];
ax = [2, 3; 4, 5; 6, 7; 8, 9];
for ii = 1:4

	% Find example in spreadsheet
	putative = neuron_time_timbre(index(iii(ii))).putative;
	acc_ideal = acc_sorted(iii(ii));
	CF = nat_data(index(iii(ii))).CF;
	MTF = nat_data(index(iii(ii))).MTF;
	fprintf('Example %d: CF = %0.0f, %s MTF\n', ii, CF, MTF)

	% Caclulate PSTH of all overlapping F0s
	table_data = getTimbreNeuronTable(nat_data, index(iii(ii)), 'Data', 20);
	for iplot = 1:32
		indices = (iplot-1)*20 + (1:20);
		mean_PSTH(iplot,:) = mean(table_data{indices,:});
	end
	basson_PSTH = mean_PSTH(1:2:32,1:15);
	oboe_PSTH = mean_PSTH(2:2:32,1:15);
	max_rate = max(max([basson_PSTH oboe_PSTH]));

	% Plot PSTHs
	h(ax(ii, 1)) = subplot(3, 4, ax(ii,1));
	plot(basson_PSTH', 'LineWidth', 1, 'Color', colors{2});
	xlim([1 15]);
	grid on;
	ylim([0 max_rate])
	xticklabels([])
	if ismember(ii, [1, 2])
		ylabel('Avg. Rate (sp/s)       ')
	end
	title(['Accuracy = ' num2str(round(acc_ideal*100, 2)) '%'])
	set(gca, 'fontsize', fontsize)
	box off
	if ii == 2
		hleg = legend('Bassoon', 'box', 'off', 'Location','northeast');
		hleg.ItemTokenSize = [6, 6];
	end

	h(ax(ii, 2)) = subplot(3, 4, ax(ii, 2));
	plot(oboe_PSTH', 'LineWidth', 1, 'Color', colors{1});
	xlim([1 15]);
	grid on;
	ylim([0 max_rate])
	xticks(0:3:15)
	xticklabels(0:60:300)
	set(gca, 'fontsize', fontsize)
	if ii == 2
		hleg = legend('Oboe', 'box', 'off', 'Location','northeast');
		hleg.ItemTokenSize = [6, 6];
	end
	xlabel('Time (ms)')
	box off
	
end


%% Scatter with rate 
accuracy_time = accuracy_all(16,:);

% Load in rate models
filepath = fullfile(base, 'model_comparisons','Neuron_Rate_Timbre_All.mat');
load(filepath, "neuron_rate_timbre")
accuracy_rate = [neuron_rate_timbre.accuracy];

h(10) = subplot(3, 4, 10);
hold on
scatter(accuracy_rate, accuracy_time, scattersize, 'filled', 'MarkerEdgeColor','k', ...
	'MarkerFaceAlpha', 0.5, 'MarkerFaceColor',colorsTimbre)
plot([0, 1], [0, 1], 'k')
plot([0 0.50], [0.50, 0.50], 'color', [0.4 0.4 0.4])
plot([0.50 0.50], [0, 0.50], 'color', [0.4 0.4 0.4])
xlim([0.40 0.90])
ylim([0.40 0.90])

% mdl = fitlm(accuracy_rate, accuracy_time);
% x = linspace(0, 1, 20);
% y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
% plot(x, y, ':k')
yticks(0:0.05:1)
xticks(0:0.05:1)
% hleg = legend('', 'Unity','Chance','', ...
% 	sprintf('y = %0.2f*x+%0.2f, p=%0.04f', mdl.Coefficients{2, 1}, ...
% 	mdl.Coefficients{1, 1},mdl.Coefficients{2,4}), 'fontsize', legsize, ...
% 	'location','best', 'box', 'off');
hleg = legend('', 'Unity','Chance','', 'fontsize', legsize, ...
	'location','best', 'box', 'off');
hleg.ItemTokenSize = [8, 8];
set(gca, 'fontsize', fontsize)
grid on
xlabel('Rate Accuracy')
ylabel('Timing Accuracy (20 ms bin)')

% Calculate and display the correlation coefficient between rate and timing accuracies
correlationCoeff = corr(accuracy_rate', accuracy_time');
fprintf('Correlation coefficient between rate and timing accuracies: %0.4f\n', correlationCoeff);

rateBetterCount = sum(accuracy_rate > accuracy_time);
timingBetterCount = sum(accuracy_time > accuracy_rate);
percentageRateBetter = (rateBetterCount / num_data) * 100;
percentageTimingBetter = (timingBetterCount / num_data) * 100;
fprintf('Rate has higher accuracy: %0.4f\n', percentageRateBetter);
fprintf('Timing has higher accuracy: %0.4f\n', percentageTimingBetter);

%[h,p,ci,stats] = ttest(accuracy_rate, accuracy_time);

%% Scatter with shuffled 

load(fullfile(base, 'model_comparisons_revised', 'Neuron_Time_Timbre_Distances_Shuffled.mat'), ...
	"neuron_time_timbre")

num_data = size(neuron_time_timbre, 2);
min_dis = [0.1 0.25 0.5 0.75 1 1.5 2 2.5 3 3.5 4 4.5 5 10 15 20 30 40 50 100 150];
accuracy_shuffled = zeros(length(min_dis), num_data);
for i = 1:num_data
	accuracy_shuffled(:,i) = [neuron_time_timbre(:,i).accuracy];
end
accuracy_rate = accuracy_shuffled(16,:);
accuracy_time = accuracy_all(16,:);

h(11) = subplot(3, 4, 11);
hold on
scatter(accuracy_rate, accuracy_time, scattersize, 'filled', 'MarkerEdgeColor','k', ...
	'MarkerFaceAlpha', 0.5, 'MarkerFaceColor',colorsTimbre)
plot([0, 1], [0, 1], 'k')
plot([0 0.50], [0.50, 0.50], 'color', [0.4 0.4 0.4])
plot([0.50 0.50], [0, 0.50], 'color', [0.4 0.4 0.4])
xlim([0.40 0.90])
ylim([0.40 0.90])

% mdl = fitlm(accuracy_rate, accuracy_time);
% x = linspace(0, 1, 20);
% y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
% plot(x, y, ':k')
yticks(0:0.05:1)
xticks(0:0.05:1)
% hleg = legend('', 'Unity','Chance','', ...
% 	sprintf('y = %0.2f*x+%0.2f, p=%0.04f', mdl.Coefficients{2, 1}, ...
% 	mdl.Coefficients{1, 1},mdl.Coefficients{2,4}), 'fontsize', legsize, ...
% 	'location','best', 'box', 'off');
hleg = legend('', 'Unity','Chance','', ...
	'fontsize', legsize, ...
	'location','best', 'box', 'off');
hleg.ItemTokenSize = [8, 8];
set(gca, 'fontsize', fontsize)
grid on
xlabel('Shuffled Timing Accuracy (20 ms bin)')
ylabel('Timing Accuracy (20 ms bin)')

% Calculate and display the correlation coefficient between rate and timing accuracies
correlationCoeff = corr(accuracy_rate', accuracy_time');
fprintf('Correlation coefficient between rate and timing accuracies: %0.4f\n', correlationCoeff);

% Calculate percentage where rate is better or timing is better
% Calculate the percentage of cases where rate accuracy is better than timing accuracy
rateBetterCount = sum(accuracy_rate > accuracy_time);
timingBetterCount = sum(accuracy_time > accuracy_rate);
percentageRateBetter = (rateBetterCount / num_data) * 100;
percentageTimingBetter = (timingBetterCount / num_data) * 100;
fprintf('Timing has higher accuracy: %0.4f\n', percentageRateBetter);
fprintf('Shuffled iming has higher accuracy: %0.4f\n', percentageTimingBetter);

[p12, ~] = ranksum(accuracy_rate, accuracy_time);  

%% Get stats 

% CF
CFs = [neuron_rate_timbre.CF];
mdl = fitlm(CFs, accuracy_time);
x = linspace(300, 14000, 50);
y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
fprintf('p=%0.04f\n',mdl.Coefficients{2,4});

%  MTF Distributions
MTFs = {neuron_rate_timbre.MTF};
MTF_types = unique(MTFs);
[weights_ordered, order_ind] = sort(accuracy_time);
hold on
accur_all_MTF = NaN(4, 123);
for iMTF = 1:4
	if iMTF == 4
		ind = strcmp(MTFs, MTF_types{iMTF}) | strcmp(MTFs, MTF_types{iMTF+1});
	else
		ind = strcmp(MTFs, MTF_types{iMTF});
	end
	[weights_ordered, order_ind] = sort(accuracy_time(ind));
	num_units = length(weights_ordered);

	mean_vals(iMTF) = mean(weights_ordered);
	std_vals(iMTF) = std(weights_ordered)/sqrt(length(weights_ordered));
	
	accur_all_MTF(iMTF, 1:num_units) = weights_ordered;

end

[p, tbl, stats] = kruskalwallis(accur_all_MTF');
[c, m, h, gnames] = multcompare(stats, 'CType', 'hsd');

%% Rearrange 

left = [0.08 0.44 0.74 0.58];
bottom = [0.1 0.55 0.8];
width = 0.24;
height = 0.07;

set(h(1), 'position', [left(1) bottom(2) 0.27 0.4])

set(h(2), 'position', [left(2) bottom(2)+height width height])
set(h(3), 'position', [left(2) bottom(2) width height])

set(h(4), 'position', [left(2) bottom(3)+height width height])
set(h(5), 'position', [left(2) bottom(3) width height])

set(h(6), 'position', [left(3) bottom(2)+height width height])
set(h(7), 'position', [left(3) bottom(2) width height])

set(h(8), 'position', [left(3) bottom(3)+height width height])
set(h(9), 'position', [left(3) bottom(3) width height])

set(h(10), 'position', [left(1) bottom(1) 0.4 0.32])
set(h(11), 'position', [left(4) bottom(1) 0.4 0.32])

%% Label plots

labelleft= [0 0.3 0.55];
labelbottom = linspace(0.28, 0.94, 2);

annotation('textbox',[0.02 0.95 0.071 0.058],...
	'String','A','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[0.35 0.95 0.071 0.058],...
	'String','B','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[0.02 0.42 0.071 0.058],...
	'String','C','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[0.52 0.42 0.071 0.058],...
	'String','D','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');


%% Save figure

if save_fig == 1
	filename = 'fig4_4_timbre_neuron_timing';
	save_figure(filename)
end



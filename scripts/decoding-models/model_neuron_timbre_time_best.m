%% model_neuron_rate_F0
clear

%% Load in data

[base, datapath, ~, ~] = getPathsNT();
load(fullfile(base,'model_comparisons_revised', 'Data_NT_3.mat'), 'nat_data')

%% Shape data into model input

[sesh, num_data] = getTimbreSessions(nat_data);
neuron_time_timbre = struct();

min_dis = [0.1 0.25 0.5 0.75 1 1.5 2 2.5 3 3.5 4 4.5 5 10 15 20 30 40 50 100 150 300];
accuracy_all = zeros(length(min_dis), num_data);

for ind = 1:num_data
	index = sesh(ind);
	parfor idist = 1:length(min_dis)
		min_dis2 = min_dis(idist);

		% Separate out training and testing data
		T = getTimbreNeuronTable(nat_data, index, 'Data', min_dis2);
		[T_train, T_test] = splitData_Reps(T);

		% Call SVM
		[trainedClassifier, validationAccuracy, validationPredictions] = ...
			trainClassifierTimeTimbre(T_train);

		% To make predictions with the returned 'trainedClassifier' on new data
		[yfit,scores] = trainedClassifier.predictFcn(T_test);
		C = confusionmat(T_test.Response, yfit);
		accuracy_all(idist, ind) = sum(diag(C))/sum(C, 'all');

		% Set up struct to save data
		neuron_time_timbre(idist, ind).putative = nat_data(index).putative;
		neuron_time_timbre(idist, ind).accuracy = accuracy_all(idist, ind);
		neuron_time_timbre(idist, ind).C = C;
		neuron_time_timbre(idist, ind).min_dis = min_dis2;


	end
	fprintf('%d/%d, %0.2f%% done!\n', ind, ...
		num_data, ind/num_data*100)
end
save(fullfile(base, 'model_comparisons_revised', 'Neuron_Time_Timbre_Distances_L1.mat'), ...
	"neuron_time_timbre", '-v7.3')
% save(fullfile(base, 'model_comparisons', 'Model_N_Time_Timbre_All.mat'), ...
% 	"neuron_time_timbre", '-v7.3')

%% Plot accuracy of each neuron

[base, ~, ~, ~] = getPathsNT();
load(fullfile(base, 'model_comparisons_revised', 'Neuron_Time_Timbre_Distances_L1.mat'), ...
	"neuron_time_timbre")

num_data = size(neuron_time_timbre, 2);
min_dis = [0.1 0.25 0.5 0.75 1 1.5 2 2.5 3 3.5 4 4.5 5 10 15 20 30 40 50 100 150 300];
accuracy_all = zeros(length(min_dis), num_data);
for i = 1:num_data
	accuracy_all(:,i) = [neuron_time_timbre(:,i).accuracy];
end

%% Get best distance for each neuron

for i = 1:num_data
	[best_acc(i), best_ind(i)] =  max(accuracy_all(:,i));
	best_dist(i) = min_dis(best_ind(i));
end

figure
histogram(best_ind)
hold on
xline(mean(best_ind), 'r')
xline(median(best_ind), '--r')
best_dist_all = min_dis(median(best_ind));
xticks(1:21)
xticklabels(min_dis)
xlim([0.5 21.5])
xlabel('PSTH Bin Width (ms)')
legend('', 'mean', 'median')
title('PSTH bin that resulted in max accuracy for each neuron')
ylabel('# Neurons')

% 10 ms is the best 
% This means that only F0s < 100 Hz can be used 

%%

figure
tiledlayout(7, 3, 'TileSpacing','tight', 'Padding','compact')
for i = 1:length(min_dis)
	nexttile
	histogram(accuracy_all(i,:)*100,51)
	mean_F0 = mean(accuracy_all(i,:));
	hold on
	xline(50, 'k')
	xline(mean_F0*100, 'r', 'LineWidth',1)
	ylabel('# Neurons')
	xlabel('Prediction Accuracy (%)')
	xlim([0 100])
	mean_all = mean(accuracy_all(i,:), 'all');
	fprintf('Mean for all = %0.4f\n', mean_all)
	title(['Min Dist = ' num2str(min_dis(i)) ', Mean = ' num2str(round(mean_all*100))])
end

%% Plot another one 

accuracy_new = accuracy_all-accuracy_all(1,:);

figure
plot(1:21, accuracy_new, 'k')
accuracy_mean = mean(accuracy_new, 2);
hold on
plot(1:21, accuracy_mean, 'r', 'linewidth', 2)
xticks(1:21)
xticklabels(min_dis)
xlim([0.5 21.5])
xlabel('PSTH Bin Width (ms)')
xline(16, 'r')
ylabel('Accuracy normalized by 0.1 ms bin accuracy')

%% Statistics 


figure
boxplot(accuracy_all')
xticks(1:21)
xticklabels(min_dis)
anova(accuracy_all')
[p,tbl,stats] = anova1(accuracy_all');
results = multcompare(stats);
yticklabels(min_dis)
ylabel('PSTH bin width (ms)')
yticklabels(fliplr(min_dis))
grid on

%% Compare with rate accuracies! 
legsize = 12;
fontsize = 12;

accuracy_time = accuracy_all(22,:);

% Load in rate models
filepath = fullfile(base, 'model_comparisons','Neuron_Rate_Timbre_All.mat');
load(filepath, "neuron_rate_timbre")
accuracy_rate = [neuron_rate_timbre.accuracy];

figure
hold on
scatter(accuracy_rate, accuracy_time, 30, 'filled', 'MarkerEdgeColor','k', ...
	'MarkerFaceAlpha', 0.5)
plot([0, 1], [0, 1], 'k')
plot([0 0.50], [0.50, 0.50], 'color', [0.4 0.4 0.4])
plot([0.50 0.50], [0, 0.50], 'color', [0.4 0.4 0.4])
xlim([0.40 0.90])
ylim([0.40 0.90])

mdl = fitlm(accuracy_rate, accuracy_time);
x = linspace(0, 1, 20);
y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
plot(x, y, ':k')
yticks(0:0.05:1)
xticks(0:0.05:1)
hleg = legend('', 'Unity','Chance','', ...
	sprintf('y = %0.2f*x+%0.2f, p=%0.04f', mdl.Coefficients{2, 1}, ...
	mdl.Coefficients{1, 1},mdl.Coefficients{2,4}), 'fontsize', legsize, ...
	'location','best', 'box', 'off');
hleg.ItemTokenSize = [8, 8];
title('Rate vs Timing (300ms bins) Comparison')
set(gca, 'fontsize', fontsize)
grid on
xlabel('Rate')
ylabel('Timing')

% Calculate and display the correlation coefficient between rate and timing accuracies
correlationCoeff = corr(accuracy_rate', accuracy_time');
fprintf('Correlation coefficient between rate and timing accuracies: %0.4f\n', correlationCoeff);

% Calculate percentage where rate is better or timing is better
% Calculate the percentage of cases where rate accuracy is better than timing accuracy
rateBetterCount = sum(accuracy_rate > accuracy_time);
timingBetterCount = sum(accuracy_time > accuracy_rate);
percentageRateBetter = (rateBetterCount / num_data) * 100;
percentageTimingBetter = (timingBetterCount / num_data) * 100;
fprintf('Rate has higher accuracy: %0.4f\n', percentageRateBetter);
fprintf('Timing has higher accuracy: %0.4f\n', percentageTimingBetter);

[h,p,ci,stats] = ttest(accuracy_rate, accuracy_time);

%% Plot PSTHs of 20 ms for high, med, and low accuracy examples 

accuracy_time = accuracy_all(16,:);
[acc_sorted, index] = sort(accuracy_time, 'descend');
iii = 192;
putative = neuron_time_timbre(index(iii)).putative;
acc_ideal = acc_sorted(iii);

% Find example in spreadsheet
s_ind = strcmp(sessions.Putative_Units, putative);
CF = sessions.CF(s_ind);

% Plot PSTH of all overlapping F0s
table_data = getTimbreNeuronTable(nat_data, index(iii), 'Data', 20);
for iplot = 1:32
	indices = (iplot-1)*20 + (1:20);
	mean_PSTH(iplot,:) = mean(table_data{indices,:});
end

% Plot PSTHs
basson_PSTH = mean_PSTH(1:2:32,1:15);
oboe_PSTH = mean_PSTH(2:2:32,1:15);
max_rate = max(max([basson_PSTH oboe_PSTH]));

figure
tiledlayout(2, 1, 'TileSpacing','none')
nexttile
plot(basson_PSTH', 'LineWidth', 1, 'Color', 'r');
xlim([1 15]);
grid on;
ylim([0 max_rate])
ylabel('Bassoon')
title(['Accuracy = ' num2str(round(acc_ideal*100, 2)) '%'])

nexttile
plot(oboe_PSTH', 'LineWidth', 1, 'Color', 'b');
xlim([1 15]);
grid on;
ylim([0 max_rate])
xticks(0:1:15)
xticklabels(0:20:300)
ylabel('Oboe')
xlabel('Time (ms)')

%% Random temporal shuffling test 



%% model_neuron_timbre_time_shuffled
clear

%% Load in data

[base, datapath, ~, ~] = getPathsNT();
load(fullfile(base,'model_comparisons_revised', 'Data_NT_3.mat'), 'nat_data')

%% Shape data into model input

[sesh, num_data] = getTimbreSessions(nat_data);
neuron_time_timbre = struct();

min_dis = 20;
accuracy_all = zeros(length(min_dis), num_data);

for ind = 1:num_data
	timerVal = tic;
	index = sesh(ind);
	parfor idist = 1:length(min_dis)
		min_dis2 = min_dis(idist);
		

		% Separate out training and testing data
		T = getTimbreNeuronTable(nat_data, index, 'Data', min_dis2);

		% Shuffle the bins in the trial (shuffle each row independently)
		data_mat2 = table2array(T(:,1:end-1));
		data_mat = zeros(640, 300/min_dis);
		for ind2 = 1:640
			row = data_mat2(ind2, :);
			shuffled_data = row(randperm(numel(row)));
			data_mat(ind2, :) = shuffled_data;
		end
		T_shuffled = array2table(data_mat);
		T_shuffled.Response = repmat([ones(20,1); ones(20, 1)*2], 16, 1);

		% Shuffle the bins in the trial (shuffle all rows together)
		data_mat2 = table2array(T(:,1:end-1));
		data_mat = zeros(640, 300/min_dis);
		for ind2 = 1:640
			row = data_mat2(ind2, :);
			shuffled_data = row(randperm(numel(row)));
			data_mat(ind2, :) = shuffled_data;
		end
		T_shuffled = array2table(data_mat);
		T_shuffled.Response = repmat([ones(20,1); ones(20, 1)*2], 16, 1);

		% Split into training and testing
		[T_train, T_test] = splitData_Reps(T_shuffled);

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
	time = toc(timerVal);
	fprintf('%d/%d, %0.2f%% done! Time=%0.2f seconds\n',...
		ind, num_data, ind/num_data*100,time)
end
save(fullfile(base, 'model_comparisons_revised', 'Neuron_Time_Timbre_Distances_Shuffled.mat'), ...
	"neuron_time_timbre", '-v7.3')
% % save(fullfile(base, 'model_comparisons', 'Model_N_Time_Timbre_All.mat'), ...
% % 	"neuron_time_timbre", '-v7.3')

%%

figure
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


%% Compare with rate accuracies! 
legsize = 12;
fontsize = 12;

accuracy_time = accuracy_all;

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
title('Rate vs Timing (20ms bins) Comparison')
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

%% Plot PSTHs of 20 ms for high, med, and low accuracy examples 

accuracy_time = accuracy_all;
[acc_sorted, index] = sort(accuracy_time, 'descend');
iii = 5;
putative = neuron_time_timbre(index(iii)).putative;
acc_ideal = acc_sorted(iii);

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

%% Load in actual results and compare! 

[base, ~, ~, ~] = getPathsNT();
load(fullfile(base, 'model_comparisons_revised', 'Neuron_Time_Timbre_Distances.mat'), ...
	"neuron_time_timbre")
legsize = 12;
fontsize = 12;

num_data = size(neuron_time_timbre, 2);
min_dis = [0.1 0.25 0.5 0.75 1 1.5 2 2.5 3 3.5 4 4.5 5 10 15 20 30 40 50 100 150];
accuracy_all = zeros(length(min_dis), num_data);
for i = 1:num_data
	accuracy_all(:,i) = [neuron_time_timbre(:,i).accuracy];
end
accuracy_rate = accuracy_all(16,:);

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
title('Timing vs Shuffled Timing (20ms bins) Comparison')
set(gca, 'fontsize', fontsize)
grid on
xlabel('Timing')
ylabel('Shuffled Timing')

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
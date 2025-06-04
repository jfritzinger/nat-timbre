%% plot_pop_rate_F0
clear

%% Load data 
%target = 'Bassoon';
target = 'Oboe';
[base, ~, ~, ~] = getPathsNT();
load(fullfile(base, 'model_comparisons', ['Pop_Rate_F0_' target '3.mat']),...
	"pop_rate_F0")
load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')

figure
tiledlayout(3, 4)
scattersize = 5;
titlesize = 9;
fontsize = 8;
labelsize = 12;
legsize = 7;


%% Plot confusion matrix 

% Compute classification using all data 
nexttile
confusionchart(pop_rate_F0.C)

% Calculate accuracy
chart = confusionchart(pop_rate_F0.response,pop_rate_F0.validationPredictions); % Generate confusion chart
confusionMatrix = chart.NormalizedValues; % Get the normalized confusion matrix
accuracy = sum(diag(confusionMatrix)) / sum(confusionMatrix(:)); % Calculate accuracy
title(sprintf('Accuracy = %0.2f%%', accuracy(1)*100))

% % Accurate % 202 / 800
% num_acc =  sum(diag(confusionMatrix));
% 
% % Accuracy within +/- 1 semitone % 13 / 800
% num_acc1 =  sum(diag(confusionMatrix, 1)) + sum(diag(confusionMatrix, -1));
% 
% % Accuracy within +/- 2 semitones % 319 / 800
% num_acc2 =  sum(diag(confusionMatrix, 2)) + sum(diag(confusionMatrix, -2));
% 
% % Accuracy within +/- 3 semitones % 319 / 800
% num_acc3 =  sum(diag(confusionMatrix, 3)) + sum(diag(confusionMatrix, -3));
% 
% % Accuracy within +/- 4 semitones % 319 / 800
% num_acc4 =  sum(diag(confusionMatrix, 4)) + sum(diag(confusionMatrix, -4));

%% Shuffled accuracy

nexttile
edges = linspace(0, 1, 41);
hold on
histogram(pop_rate_F0.shuffled_accuracy, edges)
xline(accuracy, '--r')
xline(0.5, 'k')
xlim([0 0.5])
xlabel('Accuracy')
ylabel('# Trials')
hleg = legend(['Shuffled' newline 'Data'], ['Model' newline 'Accuracy']...
	, 'Chance', 'fontsize', legsize, 'location', 'best', 'box', 'off');
hleg.ItemTokenSize = [8, 8];
set(gca, 'fontsize', fontsize)
grid on

%% Importance vs timing results 
beta_mean = pop_rate_F0.imp.ImportanceMean;
beta_std = pop_rate_F0.imp.ImportanceStandardDeviation;

load(fullfile(base, 'model_comparisons',['Neuron_Time_F0_' target '.mat']),...
	"neuron_time_F0")
accuracy_rate = [neuron_time_F0.accuracy]*100;

nexttile
hold on
scatter(accuracy_rate, abs(beta_mean),scattersize, 'filled', 'MarkerEdgeColor','k')
xlabel('Single-Unit Rate Accuracy')
ylabel('|Beta Weights|')

mdl = fitlm(accuracy_rate, abs(beta_mean));
x = linspace(0, 100, 20);
y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
plot(x, y, 'r')
hleg = legend('Neuron',sprintf('p=%0.04f',mdl.Coefficients{2,4}), ...
	'fontsize', legsize, 'location', 'northwest', 'box', 'off');
hleg.ItemTokenSize = [8, 8];
set(gca, 'fontsize', fontsize)
grid on


%% CF vs importance

CFs = [pop_rate_F0.CF];

nexttile
scatter(CFs, beta_mean, scattersize, 'filled', 'MarkerEdgeColor','k');
hold on
set(gca, 'xscale', 'log')

mdl = fitlm(CFs, beta_mean);
x = linspace(300, 14000, 50);
y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
plot(x, y, 'r')
hleg = legend('Neuron', ...
	sprintf('p=%0.04f',mdl.Coefficients{2,4}), 'fontsize', legsize, ...
	'location', 'northwest', 'box', 'off');
hleg.ItemTokenSize = [8, 8];
ylabel('Beta Weights')
xlabel('CFs')
xticks([200 500 1000 2000 5000 10000])
xticklabels([0.2 0.5 1 2 5 10])
set(gca, 'fontsize', fontsize)
grid on


%% MTF  vs importance

MTFs = pop_rate_F0.MTF;
MTF_types = unique(MTFs);
nexttile
[weights_ordered, order_ind] = sort(beta_mean);
hold on
for iMTF = 1:4
	if iMTF == 4
		ind = strcmp(MTFs, MTF_types{iMTF}) | strcmp(MTFs, MTF_types{iMTF+1});
	else
		ind = strcmp(MTFs, MTF_types{iMTF});
	end
	[weights_ordered, order_ind] = sort(beta_mean(ind));
	num_units = length(weights_ordered);

	swarmchart(ones(num_units, 1)*iMTF, weights_ordered, scattersize)

	mean_vals(iMTF) = mean(weights_ordered);
	std_vals(iMTF) = std(weights_ordered)/sqrt(length(weights_ordered));
end
errorbar(1:4, mean_vals, std_vals, 'k')
xticks(1:4)
xticklabels({'BE', 'BS', 'F', 'H'})
xlabel('MTF Groups')
ylabel('Beta Weights')

tableMTF = table(MTFs', beta_mean);
anova(tableMTF, 'beta_mean')
% [~,~,stats] = anova1(beta_mean, MTFs);
% [c,~,~,gnames] = multcompare(stats);
set(gca, 'fontsize', fontsize)
grid on

%% RVF vs importance
nexttile
if strcmp(target, 'Bassoon')
	sesh = find(~cellfun(@isempty, {nat_data.bass_rate}));
else
	sesh = find(~cellfun(@isempty, {nat_data.oboe_rate}));
end
PCA_scores = [nat_data.RVF_PC2];
PC2_score = PCA_scores(sesh)';

% Plot PCA2 score vs beta weights 
scatter(PC2_score, beta_mean, scattersize, 'filled', 'MarkerEdgeColor','k')
xlabel('PC2 RVF score')
ylabel('Beta Weights')
hold on

mdl = fitlm(PC2_score, beta_mean);
x = linspace(-2, 1, 20);
y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
plot(x, y, 'r')
hleg = legend('Neuron', ...
	sprintf('p=%0.04f',mdl.Coefficients{2,4}), 'fontsize', legsize, ...
	'location', 'northwest', 'box', 'off');
hleg.ItemTokenSize = [8, 8];
set(gca, 'fontsize', fontsize)
grid on

%% Subset CF groups 
% 
% load(fullfile(base, 'model_comparisons', 'Pop_Rate_Timbre_Subset_CF.mat'), ...
% 	"accuracies", "CF_names")
% 
% hold on
% for iCF = 1:7
% 	swarmchart(ones(size(accuracies, 2), 1)*iCF, accuracies(iCF,:), scattersize)
% end
% xlabel('CF Groups')
% xticks(1:7)
% xticklabels(CF_names)
% max_acc = max(accuracies, [], 2);
% plot(1:7, max_acc, 'k')
% mean_acc = median(accuracies, 2);
% plot(1:7, mean_acc, 'k')
% set(gca, 'fontsize', fontsize)
% ylabel('Model Accuracy')
% grid on
% 
% %% I. Subset MTF groups 
% 
% load(fullfile(base, 'model_comparisons', 'Pop_Rate_Timbre_Subset_MTF.mat'), ...
% 	"accuracy_all", "MTF_names")
% 
% hold on
% for iCF = 1:5
% 	swarmchart(ones(size(accuracy_all, 2), 1)*iCF, accuracy_all(iCF,:), scattersize)
% end
% xlabel('MTF Groups')
% xticks(1:5)
% xticklabels(MTF_names)
% ylim([0.78 1])
% max_acc = max(accuracy_all, [], 2);
% plot(1:5, max_acc, 'k')
% mean_acc = median(accuracy_all, 2);
% plot(1:5, mean_acc, 'k')
% set(gca, 'fontsize', fontsize)
% ylabel('Model Accuracy')
% grid on
% 
% %% G. Subset accuracies 
% 
% load(fullfile(base, 'model_comparisons', 'Pop_Rate_Timbre_Best.mat'), ...
% 	"pop_rate_timbre_best")
% 
% nneurons = [1:4 5:5:50];
% accuracy_bad = [pop_rate_timbre_best(15:28).accuracy];
% plot(nneurons, accuracy_bad);
% 
% hold on 
% accuracy_good = [pop_rate_timbre_best(1:14).accuracy];
% plot(nneurons, accuracy_good);
% xlabel('# Neurons in Model')
% ylabel('Model Accuracy')
% grid on
% box off
% hleg = legend('Worst', 'Best', 'fontsize', legsize, 'location', 'best', 'box', 'off');
% hleg.ItemTokenSize = [8, 8];
% set(gca, 'fontsize', fontsize)

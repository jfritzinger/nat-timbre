%% Plot model evaluations 
clear

%% Load in model

[base, datapath, savepath, ppi] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Pop_Rate_Timbre_All.mat'), ...
	"pop_rate_timbre")

%% Plot histogram of accuracy vs shuffled
figure
tiledlayout(1, 2)

% Compute classification using all data 
nexttile
C = confusionmat(pop_rate_timbre.response, pop_rate_timbre.predictions);
confusionchart(C)
title(sprintf('Accuracy = %0.2f%%', pop_rate_timbre.accuracy*100))

nexttile
edges = linspace(0, 100, 101);
hold on
histogram(pop_rate_timbre.shuffled_accuracy*100)
xline(pop_rate_timbre.accuracy*100, '--r')
xline(50, 'k')
xlim([0 101])
xlabel('Accuracy (%)')
ylabel('# Trials')
title('Accuracy for Shuffled Data')
legend(...
	sprintf('Mean = %0.0f%%', mean(pop_rate_timbre.shuffled_accuracy*100)))

%% Plots example rates for example F0s

% % Plot all rates/VS 
% figure
% tiledlayout(2, 1)
% mean_bass = mean(data_mat(1:20,:),1);
% mean_oboe = mean(data_mat(21:40,:),1);
% for ii = 1:2
% 	nexttile
% 	if ii == 1
% 		scatter(CFs/1000, mean_bass, 'blue')
% 		title('Bassoon')
% 	else
% 		scatter(CFs/1000, mean_oboe, 'red')
% 		title('Oboe')
% 	end
% 	set(gca, 'xscale', 'log')
% 	xticks([0.1 0.2 0.5 1 2 5 10])
% 	xlim([0.1 10])
% 	ylabel('Norm Rate')
% 	xlabel('CF (kHz)')
% 	grid on
% end

%% Initial beta weight investigation

beta_weights = pop_rate_timbre.trainedClassifier.ClassificationSVM.Beta;
CFs = pop_rate_timbre.CFs;

figure
tiledlayout(1, 3)

nexttile
plot(1:180, beta_weights)
hold on 
yline(0)
xlabel('Neuron #')
ylabel('Beta Weight')

nexttile
scatter(CFs/1000, beta_weights)
hold on 
yline(0)
xticks([0.1 0.2 0.5 1 2 5 10])
xlim([0.2 10])
xlabel('CF')
ylabel('Beta Weight')
set(gca, 'xscale', 'log')

nexttile
scatter(CFs/1000, abs(beta_weights))
hold on 
yline(0)
xticks([0.1 0.2 0.5 1 2 5 10])
xlim([0.2 10])
xlabel('CF')
ylabel('abs(Beta Weight)')
set(gca, 'xscale', 'log')

%% Get top/bottom 3 neurons with highest abs(beta_weights)

[~,originalpos] = sort(abs(beta_weights), 'descend' );
best_ind=originalpos(1:3);
putatives_best = pop_rate_timbre.putative(best_ind);
best_ind_15 = originalpos(1:15);

[~,originalpos] = sort(abs(beta_weights), 'ascend' );
worst_ind=originalpos(1:3);
putatives_worst = pop_rate_timbre.putative(worst_ind);
worst_ind_15 = originalpos(1:15);

%% Search for patterns in MTF in beta weights 

beta_weights = pop_rate_timbre.trainedClassifier.ClassificationSVM.Beta;
CFs = pop_rate_timbre.CFs;
MTFs = pop_rate_timbre.MTF;
MTF_types = unique(MTFs);

[weights_ordered, order_ind] = sort(beta_weights);
figure
bar(weights_ordered)

figure
for iMTF = 1:5

	nexttile
	ind = strcmp(MTFs, MTF_types{iMTF});
	[weights_ordered, order_ind] = sort(beta_weights(ind));
	bar(weights_ordered)
	title(MTF_types{iMTF})
	ylim([-2 2])
end

figure
hold on
for iMTF = 1:5
	ind = strcmp(MTFs, MTF_types{iMTF});
	[weights_ordered, order_ind] = sort(beta_weights(ind));
	num_units = length(weights_ordered);

	swarmchart(ones(num_units, 1)*iMTF, weights_ordered)
	ylim([-2 2])

	mean_vals(iMTF) = mean(weights_ordered);
	std_vals(iMTF) = std(weights_ordered)/sqrt(length(weights_ordered));
end
errorbar(1:5, mean_vals, std_vals, 'k')
xticks(1:5)
xticklabels(MTF_types)
ylabel('Beta Weights')

tableMTF = table(MTFs', beta_weights);
anova(tableMTF, 'beta_weights')
%[~,~,stats] = anova1(beta_weights, MTFs);
%[c,~,~,gnames] = multcompare(stats);


%% CFs 

figure
CF_groups = [0, 2000; 2000, 4000; 4000, 14000];
CF_names = {'Low', 'Medium', 'High'};
for iCF = 1:3

	nexttile
	ind = CFs > CF_groups(iCF, 1) & CFs < CF_groups(iCF, 2);
	[weights_ordered, order_ind] = sort(beta_weights(ind));
	bar(weights_ordered)
	title(CF_names{iCF})
	ylim([-2 2])
end

figure
hold on
mean_vals = zeros(1, 3);
std_vals = zeros(1,3);
for iCF = 1:3
	ind = CFs > CF_groups(iCF, 1) & CFs < CF_groups(iCF, 2);
	[weights_ordered, order_ind] = sort(beta_weights(ind));
	num_units = length(weights_ordered);

	swarmchart(ones(num_units, 1)*iCF, weights_ordered)
	ylim([-2 2])

	mean_vals(iCF) = mean(weights_ordered);
	std_vals(iCF) = std(weights_ordered)/sqrt(length(weights_ordered));
end
errorbar(1:3, mean_vals, std_vals, 'k')
xticks(1:3)
xticklabels(CF_names)
ylabel('Beta Weights')

tableMTF = table(CFs', beta_weights);
anova(tableMTF, 'beta_weights')
%[~,~,stats] = anova1(beta_weights, MTFs);
%[c,~,~,gnames] = multcompare(stats);

%% Throw into a mixed-effects model and see if anything stands out 




%%
%%
% 
% figure
% tiledlayout(2, 1)
% mean_bass = mean(data_mat(1:20,:),1);
% mean_oboe = mean(data_mat(21:40,:),1);
% for ii = 1:2
% 	nexttile
% 	if ii == 1
% 		scatter(CFs/1000, mean_bass, [], beta_weights, 'filled', ...
% 			'MarkerEdgeColor','k')
% 		title('Bassoon')
% 	else
% 		scatter(CFs/1000, mean_oboe, [], beta_weights, 'filled', ...
% 			'MarkerEdgeColor','k')
% 		title('Oboe')
% 	end
% 	set(gca, 'xscale', 'log')
% 	xticks([0.1 0.2 0.5 1 2 5 10])
% 	xlim([0.2 10])
% 	ylabel('Norm Rate')
% 	xlabel('CF (kHz)')
% 	grid on
% 	colorbar
% end
% 
% figure
% tiledlayout(1, 2)
% mean_bass = mean(data_mat(1:20,:),1);
% mean_oboe = mean(data_mat(21:40,:),1);
% index = abs(beta_weights)>0.4;
% for ii = 1:2
% 	nexttile
% 	if ii == 1
% 		scatter(CFs(index)/1000, mean_bass(index), [], beta_weights(index), 'filled', ...
% 			'MarkerEdgeColor','k')
% 		title('Bassoon')
% 	else
% 		scatter(CFs(index)/1000, mean_oboe(index), [], beta_weights(index), 'filled', ...
% 			'MarkerEdgeColor','k')
% 		title('Oboe')
% 	end
% 	set(gca, 'xscale', 'log')
% 	xticks([0.1 0.2 0.5 1 2 5 10])
% 	xlim([0.2 10])
% 	ylim([0 80])
% 	ylabel('Norm Rate')
% 	xlabel('CF (kHz)')
% 	grid on
% 	colorbar
% end

%%
% 
% figure
% tiledlayout(1, 2)
% nexttile
% scatter(mean_bass, mean_oboe, [], beta_weights, 'filled','MarkerEdgeColor','k')
% hold on
% plot([0 100], [0 100], 'k')
% xlabel('Bassoon Norm Rate')
% ylabel('Oboe Norm Rate')
% xlim([0 100])
% ylim([0 100])
% colorbar
% 
% nexttile
% scatter(mean_bass(index), mean_oboe(index), 'filled','MarkerEdgeColor','k')
% hold on
% plot([0 100], [0 100], 'k')
% xlabel('Bassoon Norm Rate')
% ylabel('Oboe Norm Rate')
% xlim([0 100])
% ylim([0 100])

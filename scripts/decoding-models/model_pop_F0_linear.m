clear

%% Load in data
%target = 'Bassoon';
%target = 'Oboe';
target = 'Invariant';

base = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')
tuning = readtable(fullfile(base, 'Tuning.xlsx')); % Load in tuning

%% Get data into proper matrix

F0s = getF0s(target);
F0s = log10(F0s);
[sesh, num_data] = getF0Sessions(nat_data, target);
T = getF0PopTable(nat_data, target, sesh, F0s, num_data, 'linear');
CFs = [nat_data(sesh).CF];
MTFs = {nat_data(sesh).MTF};
putative = {nat_data(sesh).putative};


%% Fit using fitrlinear

Mdl = fitrlinear(T, 'Response','BetaTolerance',0.0001, ...
	'Learner','leastsquares', 'Lambda','auto', 'Solver','lbfgs', ...
	'KFold',5, 'CrossVal','on', 'Regularization','ridge');
pred_F0 = kfoldPredict(Mdl);

save(fullfile(base, 'model_comparisons', ['Pop_Rate_F0_Linear_' target '.mat']),...
	"pred_F0", "T", "r2", "C", "accuracy", "F0s", "Mdl")

%% Plot results
scattersize = 20;

r = corrcoef(pred_F0, T.Response);
r2 = r(1, 2)^2;
actual = T.Response; % Actual response values
mse = mean((actual - pred_F0).^2);   % Mean Squared Error
rmse = sqrt(mse);                     % Root Mean Squared Error
mae = mean(abs(actual - pred_F0));    % Mean Absolute Error

figure('Position',[31,910,910,413])
tiledlayout(2, 4)
nexttile
scatter(10.^(T.Response), 10.^(pred_F0), scattersize, 'filled', ...
		'MarkerFaceAlpha',0.5, 'MarkerFaceColor','k')
%scatter(T.Response, pred_F0, 'filled', 'MarkerFaceAlpha',0.5)
set(gca, 'xscale', 'log', 'yscale', 'log')
% ylim([57 588])
% xlim([57 588])
% xticks([100 200 500])
hold on
plot(10.^T.Response, 10.^T.Response, 'k')
%plot(T.Response, T.Response)
title(['R^2 = ' num2str(r2)])
xlabel('Actual F0s (Log10(F0s))')
ylabel('Predicted F0s (Log10(F0s))')
grid on

closest_cat = [];
for ii = 1:length(pred_F0)
	differences = abs(F0s - pred_F0(ii));
	[~, closest_column_index] = min(differences);
	closest_cat(ii) = F0s(closest_column_index);
end

% Plot confusion matrix
nexttile
C = confusionmat(T.Response, closest_cat);
confusionchart(C)

% Calculate accuracy
chart = confusionchart(T.Response,closest_cat); % Generate confusion chart
confusionMatrix = chart.NormalizedValues; % Get the normalized confusion matrix
accuracy = sum(diag(confusionMatrix)) / sum(confusionMatrix(:)); % Calculate accuracy
title(sprintf('Accuracy = %0.2f%%', accuracy*100))

% %% Residuals analysis
% figure
% residuals = actual - pred_F0;
% scatter(pred_F0, residuals);
% xlabel('Predicted Values'); 
% ylabel('Residuals');

%% Beta weight analysis

trainedModels = Mdl.Trained;  % Cell array of models (one per fold)

% Preallocate matrix for coefficients
numFeatures = size(trainedModels{1}.Beta, 1);
allBetas = zeros(numFeatures, Mdl.KFold);

% Populate coefficients from each fold
for fold = 1:Mdl.KFold
    allBetas(:, fold) = trainedModels{fold}.Beta;
end
avgBeta = mean(allBetas, 2);  % Average across folds

[~,originalpos] = sort(avgBeta, 'descend' );
ordered_mean = avgBeta(originalpos);

nexttile
bar(ordered_mean)
hold on
plot(1:length(avgBeta), ordered_mean, 'color', 'k', 'LineStyle', 'none');
%xlim([0 25])
box off
xlabel('Neuron #')
ylabel('Importance')
grid on
title('Beta Weights')

%% Analysis of beta weights and neuron qualities 
legsize = 8;

% Beta weights, CF group
nexttile
scatter(CFs, avgBeta, scattersize, 'filled', 'MarkerEdgeColor','k');
hold on
set(gca, 'xscale', 'log')

mdl = fitlm(CFs, avgBeta);
x = linspace(300, 14000, 50);
y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
plot(x, y, ':k')
hleg = legend('Neuron', ...
	sprintf('p=%0.04f',mdl.Coefficients{2,4}), 'fontsize', legsize, ...
	'location', 'northwest', 'box', 'off');
hleg.ItemTokenSize = [8, 8];
ylabel('Beta Weights')
xlabel('CFs')
xticks([200 500 1000 2000 5000 10000])
xticklabels([0.2 0.5 1 2 5 10])
grid on

%% D. Beta weights, MTF types 

nexttile
colorsMTF = {'#648FFF', '#DC267F', '#785EF0', '#FFB000'};
MTF_types = unique(MTFs);
[weights_ordered, order_ind] = sort(avgBeta);
hold on
for iMTF = 1:4
	if iMTF == 4
		ind = strcmp(MTFs, MTF_types{iMTF}) | strcmp(MTFs, MTF_types{iMTF+1});
	else
		ind = strcmp(MTFs, MTF_types{iMTF});
	end
	[weights_ordered, order_ind] = sort(avgBeta(ind));
	num_units = length(weights_ordered);

	RGB = hex2rgb(colorsMTF{iMTF});
	swarmchart(ones(num_units, 1)*iMTF, weights_ordered, scattersize, RGB)

	mean_vals(iMTF) = mean(weights_ordered);
	std_vals(iMTF) = std(weights_ordered)/sqrt(length(weights_ordered));
end
errorbar(1:4, mean_vals, std_vals, 'k')
xticks(1:4)
xticklabels({'BE', 'BS', 'F', 'H'})
xlabel('MTF Groups')
ylabel('Beta Weights')

% tableMTF = table(MTFs', avgBeta);
% anova(tableMTF, 'avgBeta')
% [~,~,stats] = anova1(avgBeta, MTFs);
% [c,~,~,gnames] = multcompare(stats);

% Beta weights vs PCA for RVF 
nexttile
PCA_scores = [nat_data.RVF_PC2];
PC2_score = PCA_scores(sesh)';

% Plot PCA2 score vs beta weights 
scatter(PC2_score, avgBeta, scattersize, 'filled', 'MarkerEdgeColor','k')
xlabel('PC2 RVF score')
ylabel('Beta Weights')
hold on

mdl = fitlm(PC2_score, avgBeta);
x = linspace(-2, 1, 20);
y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
plot(x, y, ':k')
hleg = legend('Neuron', ...
	sprintf('p=%0.04f',mdl.Coefficients{2,4}), 'fontsize', legsize, ...
	'location', 'northwest', 'box', 'off');
hleg.ItemTokenSize = [8, 8];
grid on

%% Plot rate responses for 'best' neurons 

% [~,best_ind] = sort(abs(avgBeta), 'descend');
% best_put = putative(best_ind(1:16));
% F0s_oboe = getF0s('Oboe');
% F0s_bass = getF0s('Bassoon');
% F0s_both = getF0s('Invariant');
% 
% figure
% tiledlayout(4, 4)
% colorsData = {"#0072BD", "#D95319"};
% linewidth = 1;
% 
% for ii = 1:16
% 
% 	nexttile
% 	% Get rates
% 	index = find(strcmp(best_put{ii}, {nat_data.putative}));
% 	% oboe_rate = nat_data(index).oboe_rate(ind_o);
% 	% oboe_std = nat_data(index).oboe_rate_std(ind_o);
% 	% bass_rate = nat_data(index).bass_rate(ind_b);
% 	% bass_std = nat_data(index).bass_rate_std(ind_b);
% 	oboe_rate = nat_data(index).oboe_rate;
% 	oboe_std = nat_data(index).oboe_rate_std;
% 	bass_rate = nat_data(index).bass_rate;
% 	bass_std = nat_data(index).bass_rate_std;
% 	CF = nat_data(index).CF;
% 	MTF = nat_data(index).MTF;
% 
% 	% errorbar(pitches, oboe_rate, oboe_std/sqrt(20), 'color',colorsData{1}, ...
% 	% 	'LineWidth', linewidth, 'CapSize',3);
% 	% fill([F0s_both(1) F0s_both(end) F0s_both(end) ...
% 	% 	F0s_both(1)], [0 0 max([oboe_rate; bass_rate]*1.3) ...
% 	% 	max([oboe_rate; bass_rate]*1.3)], 'k', 'FaceAlpha',0.1, ...
% 	% 	'EdgeColor','none')
% 	hold on
% 	if strcmp(target, 'Oboe')
% 		errorbar(F0s_oboe, oboe_rate, oboe_std/sqrt(20), 'color',colorsData{1}, ...
% 			'LineWidth', linewidth, 'CapSize',3);
% 	elseif strcmp(target, 'Bassoon')
% 		errorbar(F0s_bass, bass_rate, bass_std/sqrt(20), 'color',colorsData{2}, ...
% 			'LineWidth', linewidth, 'CapSize',3);
% 	else
% 		errorbar(F0s_oboe, oboe_rate, oboe_std/sqrt(20), 'color',colorsData{1}, ...
% 			'LineWidth', linewidth, 'CapSize',3);
% 		errorbar(F0s_bass, bass_rate, bass_std/sqrt(20), 'color',colorsData{2}, ...
% 			'LineWidth', linewidth, 'CapSize',3);
% 	end
% 	box off
% 	grid on
% 	xticks([50 100 250 500 1000])
% 	ylim([0 max([oboe_rate; bass_rate]*1.3)])
% 	yticks([0 25 50 75 100])
% 
% 	if ismember(ii, [1 5 9 13]) 
% 		ylabel('Avg Rate (sp/s)')
% 	end
% 
% 	if ismember(ii, 13:16)
% 		xlabel('F0 (kHz)')
% 		xticklabels([50 100 250 500 1000]/1000)
% 	else
% 		xticklabels([])
% 	end
% 
% 	set(gca, 'XScale', 'log')
% 
% end

%% Run subsets 

% Get best units
beta_weights = avgBeta;
[~,best_ind] = sort(abs(beta_weights), 'descend' );

% Find all rows with bassoon and oboe
[sesh_all, num_data_all] = getF0Sessions(nat_data, target);
CFs_all = [nat_data(sesh_all).CF];
putative = {nat_data(sesh_all).putative};
MTFs = {nat_data(sesh_all).MTF};

% Get a subset of the data based on MTF type
MTF_names = {'All', 'BE', 'BS', 'H', 'F'};
ind_all(1,:) = ones(1,num_data_all);
ind_all(2,:) = strcmp('BE', MTFs);
ind_all(3,:) = strcmp('BS', MTFs);
ind_all(4,:) = contains(MTFs,'H');
ind_all(5,:) = strcmp('F', MTFs);
ind_all = logical(ind_all);

% Get a subset of the data based on CF grouping
CF_groups = [0, 14000; 0, 2000; 2000 4000; 4000 8000];
num_subset = [1 5:5:90 100:10:200];

%% Model CF Groups 

iii = 0;
nreps = 1;
accuracy_all = NaN(length(CF_groups), length(num_subset), nreps);
for iCF = 1:length(CF_groups)
	ind = CFs_all > CF_groups(iCF, 1) & CFs_all < CF_groups(iCF, 2);
	beta_weights = avgBeta;
	[~,best_ind] = sort(abs(beta_weights(ind)), 'descend');
	for inum = 1:length(num_subset)
		timerVal = tic;
		for irep = 1:nreps
			if num_subset(inum)>length(best_ind)
				accuracy_all(iCF, inum, irep) = NaN;
			else
				% Get subset
				rand_ind = best_ind(1:num_subset(inum));
				CFs = CFs_all(rand_ind);
				sesh = sesh_all(rand_ind);
				num_data = numel(sesh);

				% Create table for model
				T = getF0PopTable(nat_data, target, sesh, F0s, num_data, 'linear');

				% Run model with kfold validation
				Mdl = fitrlinear(T, 'Response','BetaTolerance',0.0001, ...
					'Learner','leastsquares', 'Lambda','auto', 'Solver','lbfgs', ...
					'KFold',5, 'CrossVal','on', 'Regularization','ridge');
				pred_F0 = kfoldPredict(Mdl);
				r = corrcoef(pred_F0, T.Response);
				r2 = r(1, 2)^2;
				accuracy_all(iCF, inum, irep) = r2;
			end
		end
		iii = iii+1;
		timer = toc(timerVal);
		fprintf('Models took %0.2g seconds, %d/%d, %0.0f%% done\n', ...
			timer, iii, length(CF_groups)*length(num_subset), ...
			iii/(length(CF_groups)*length(num_subset))*100)
	end
end

% Plot
accuracy = mean(accuracy_all, 3);
nexttile
hold on
plot(num_subset, accuracy', 'linewidth', 2)
ylabel('Accuracy')
xlabel('# Neurons in Model')
%legend('All', 'CF = 0-2 kHz', 'CF = 2-4 kHz', 'CF = 4-8 kHz')
grid on 
box off

% Save accuracies
save(fullfile(base, 'model_comparisons', ['Pop_Rate_F0_' target '_CF.mat']), ...
	"accuracy", "CF_groups", "num_subset")

%% Model MTF Groups

iii = 0;
accuracy_all = NaN(length(MTF_names), length(num_subset), nreps);
for iMTF = 1:length(MTF_names)

	% Get index 
	ind = ind_all(iMTF,:);
	beta_weights = avgBeta;
	[~,best_ind] = sort(abs(beta_weights(ind)), 'descend');
	for inum = 1:length(num_subset)
		timerVal = tic;

		for irep = 1:nreps
			if num_subset(inum)>length(best_ind)
				accuracy_all(iMTF, inum, irep) = NaN;
			else
				% Get subset
				rand_ind = best_ind(1:num_subset(inum));
				CFs = CFs_all(rand_ind);
				sesh = sesh_all(rand_ind);
				num_data = numel(sesh);

				% Create table for model
				T = getF0PopTable(nat_data, target, sesh, F0s, num_data, 'linear');

				% Run model with kfold validation
				Mdl = fitrlinear(T, 'Response','BetaTolerance',0.0001, ...
					'Learner','leastsquares', 'Lambda','auto', 'Solver','lbfgs', ...
					'KFold',5, 'CrossVal','on', 'Regularization','ridge');
				pred_F0 = kfoldPredict(Mdl);
				r = corrcoef(pred_F0, T.Response);
				r2 = r(1, 2)^2;
				accuracy_all(iMTF, inum, irep) = r2;

			end
		end
		iii = iii+1;
		timer = toc(timerVal);
		fprintf('Models took %0.2g seconds, %d/%d, %0.0f%% done\n', ...
			timer, iii, length(MTF_names)*length(num_subset), ...
			iii/(length(MTF_names)*length(num_subset))*100)
	end
end

% Plot
accuracy = mean(accuracy_all, 3);
std_acc = std(accuracy_all, [],3);
colorsMTF = {'#1b9e77', '#648FFF', '#DC267F', '#785EF0', '#FFB000'};

nexttile
hold on
for iMTF = 1:5
	plot(num_subset, accuracy(iMTF,:), 'Color', colorsMTF{iMTF}, 'linewidth', 2)
end
ylabel('Accuracy')
xlabel('# Neurons in Model')
%legend(MTF_names)
grid on 
box off

% Save accuracies
save(fullfile(base, 'model_comparisons', ['Pop_Rate_F0_' target '_MTF.mat']), ...
	"accuracy", "MTF_names", "num_subset")
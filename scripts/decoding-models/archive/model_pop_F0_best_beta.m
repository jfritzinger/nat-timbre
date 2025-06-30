%% model_pop_VS_F0
clear

%% Load in data
%target = 'Bassoon';
%target = 'Oboe';
target = 'Invariant';

[base, ~, ~, ~] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')
load(fullfile(base, 'model_comparisons', ['Pop_Rate_F0_' target '.mat']), ...
	"pop_rate_F0")

%% Get correct output of model

% Get stimulus
F0s = getF0s(target);

% Find all rows with bassoon and oboe
[sesh_all, num_data_all] = getF0Sessions(nat_data, target);
CFs_all = [nat_data(sesh_all).CF];
putative = {nat_data(sesh_all).putative};
MTFs = {nat_data(sesh_all).MTF};

% Get indices of interest
beta_weights = pop_rate_F0.imp.ImportanceMean;
[~,originalpos] = sort(abs(beta_weights), 'descend' );
best_ind = originalpos(1:200);
[~,originalpos] = sort(abs(beta_weights), 'ascend' );
worst_ind = originalpos(1:200);

num_neurons = [1:4 5:5:50 60:10:100 150 200 1:4 5:5:50 60:10:100 150 200];
nmodels = length(num_neurons);
%poolobj = parpool();
for imodel = 1:nmodels
	timerVal = tic;

	for irep = 1:10
		if imodel < nmodels/2+1 % 6 good models
			index = best_ind;
		else % 6 bad models
			index = worst_ind;
		end
		num_index = 1:num_neurons(imodel);
		num_data = num_neurons(imodel);
		sesh = sesh_all(index(num_index));

		% Create table for model
		T = getF0PopTable(nat_data, target, sesh, F0s, num_data);

		% Train model
		[trainedClassifier1, validationAccuracy1, validationPredictions1] = ...
			trainClassifierPopRateF0(T, target);
		C = confusionmat(T.Response, validationPredictions1);
		accuracies(imodel, irep) = sum(diag(C))/sum(C, 'all');
	end
	fprintf('Models took %0.2g minutes\n', toc(timerVal)/60)
	fprintf('%d/%d, %0.2f%% done!\n', imodel, nmodels, imodel/nmodels*100)

end
%delete(poolobj);

%% Plot results

mean_acc = mean(accuracies,2);
std_acc = std(accuracies, [], 2);

figure
errorbar(num_neurons(1:nmodels/2), mean_acc(1:nmodels/2), std_acc(1:nmodels/2)/sqrt(10));
hold on
errorbar(num_neurons(nmodels/2+1:end), mean_acc(nmodels/2+1:end), std_acc(nmodels/2+1:end)/sqrt(10))
xlabel('Number of Neurons in Model')
ylabel('Accuracy')
legend('Best', 'Worst')

%% Plot outputs

save(fullfile(base, 'model_comparisons', ['Pop_Rate_F0_Best_' target '.mat']), ...
	"accuracies", "num_neurons")


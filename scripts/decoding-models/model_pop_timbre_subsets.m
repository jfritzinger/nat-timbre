%% model_pop_timbre_subsets
clear

%% Load in data

[base, ~, ~, ~] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')
load(fullfile(base, 'model_comparisons', 'Pop_Rate_Timbre_All.mat'), ...
	"pop_rate_timbre")
% load(fullfile(base, 'model_comparisons', 'Model_Pop_Rate_Timbre_All.mat'), ...
% 	"pop_rate_timbre")
% load(fullfile(base, 'model_comparisons',  'Model_NT.mat'), 'nat_model')
% nat_data = nat_model;

%% Get data into proper matrix

% Get best units
beta_weights = pop_rate_timbre.trainedClassifier.ClassificationSVM.Beta;
[~,best_ind] = sort(abs(beta_weights), 'descend' );
[~,worst_ind] = sort(abs(beta_weights), 'ascend' );

% Find all rows with bassoon and oboe
[sesh_all, num_data_all] = getTimbreSessions(nat_data);
CFs_all = [nat_data(sesh_all).CF];

%% CF group subset

% Get a subset of the data based on CF grouping
CF_groups = [0, 14000; 0, 2000; 2000 4000; 4000 8000];
num_subset = [1, 2, 3, 4, 5:5:90];
iii = 0;
accuracy_all = NaN(length(CF_groups), length(num_subset), 1);
for iCF = 1:length(CF_groups)
	ind = CFs_all > CF_groups(iCF, 1) & CFs_all < CF_groups(iCF, 2);

	beta_weights = pop_rate_timbre.trainedClassifier.ClassificationSVM.Beta;
	[~,best_ind] = sort(abs(beta_weights(ind)), 'descend' );
	for inum = 1:length(num_subset)
		timerVal = tic;

		for irep = 1
			if num_subset(inum)>length(best_ind)
				accuracy_all(iCF, inum, irep) = NaN;
			else

				% Subset!
				rand_ind = best_ind(1:num_subset(inum));
				CFs = CFs_all(rand_ind);
				sesh = sesh_all(rand_ind);
				num_data = numel(sesh);

				% Get data into table 
				T = getTimbrePopTable(nat_data, 'Rate', sesh, num_data);

				% Separate out training and testing data
				[T_train, T_test] = splitData(T); 

				% Run model with kfold validation
				[trainedClassifier, ~, ~] = trainClassifierPopRateTimbre(T_train);

				% To make predictions with the returned 'trainedClassifier' on new data
				[yfit,scores] = trainedClassifier.predictFcn(T_test);
				C = confusionmat(T_test.Response, yfit);
				accuracy = sum(diag(C))/sum(C, 'all');

				% pop_rate_timbre.trainedClassifier = trainedClassifier;
				% pop_rate_timbre.accuracy = accuracy;
				% pop_rate_timbre.predictions = predictions;
				% pop_rate_timbre.response = T.Instrument;
				% pop_rate_timbre.T = T;
				% pop_rate_timbre.CFs = CFs;
				% pop_rate_timbre.putative = putative;
				% pop_rate_timbre.sesh = sesh;
				% pop_rate_timbre.MTF = {nat_data(sesh).MTF};
				% pop_rate_timbre.oboe_rate = [nat_data(sesh).oboe_rate];
				% pop_rate_timbre.oboe_rate_std = [nat_data(sesh).oboe_rate_std];
				% pop_rate_timbre.bass_rate = [nat_data(sesh).bass_rate];
				% pop_rate_timbre.bass_rate_std = [nat_data(sesh).bass_rate_std];

				accuracy_all(iCF, inum, irep) = accuracy;
			end
		end
		iii = iii+1;
		timer = toc(timerVal);
		fprintf('Models took %0.2g seconds, %d/%d, %0.0f%% done\n', ...
			timer, iii, 4*length(num_subset), iii/(4*length(num_subset))*100)
	end
end

% Plot results 
accuracy = mean(accuracy_all, 3);
std_acc = std(accuracy_all, [],3);
figure
hold on
errorbar(num_subset, accuracy', std_acc'/sqrt(5), 'linewidth', 2)
ylabel('Accuracy')
xlabel('# Neurons in Model')
legend('All', 'CF = 0-2 kHz', 'CF = 2-4 kHz', 'CF = 4-8 kHz')
grid on 
box off

% Save accuracies
% save(fullfile(base, 'model_comparisons', 'Pop_Rate_Timbre_Subset_CF.mat'), ...
% 	"accuracy", "std_acc", "CF_groups")
save(fullfile(base, 'model_comparisons', 'Model_Pop_Rate_Timbre_Subset_CF.mat'), ...
	"accuracy", "std_acc", "CF_groups")

%% MTF Subset

% Get a subset of the data based on CF grouping
MTF_names = {'All', 'BE', 'BS', 'H', 'F'};

% Get a subset of the data based on MTF type
MTFs = {nat_data(sesh_all).MTF};
ind_all(1,:) = ones(1,num_data_all);
ind_all(2,:) = strcmp('BE', MTFs);
ind_all(3,:) = strcmp('BS', MTFs);
ind_all(4,:) = contains(MTFs,'H');
ind_all(5,:) = strcmp('F', MTFs);
ind_all = logical(ind_all);

num_subset = [1, 2, 3, 4, 5:5:90];
iii = 0;
accuracy_all = NaN(length(MTF_names), length(num_subset), 1);
for iMTF = 1:length(MTF_names)
	ind = ind_all(iMTF,:);

	beta_weights = pop_rate_timbre.trainedClassifier.ClassificationSVM.Beta;
	[~,best_ind] = sort(abs(beta_weights(ind)), 'descend' );
	for inum = 1:length(num_subset)
		timerVal = tic;

		for irep = 1
			if num_subset(inum)>length(best_ind)
				accuracy_all(iMTF, inum, irep) = NaN;
			else

				% Subset!
				rand_ind = best_ind(1:num_subset(inum));
				CFs = CFs_all(rand_ind);
				sesh = sesh_all(rand_ind);
				num_data = numel(sesh);

				% Get data into table 
				T = getTimbrePopTable(nat_data, 'Rate', sesh, num_data);

				% Run model with kfold validation
				[trainedClassifier, accuracy, predictions] = trainClassifierPopRateTimbre(T);
				accuracy_all(iMTF, inum, irep) = accuracy;
			end
		end
		iii = iii+1;
		timer = toc(timerVal);
		fprintf('Models took %0.2g seconds, %d/%d, %0.0f%% done\n', timer, ...
			iii, 5*length(num_subset), iii/(5*length(num_subset))*100)
	end
end

% Plot results 
accuracy = mean(accuracy_all, 3);
std_acc = std(accuracy_all, [],3);
colorsMTF = {'#1b9e77', '#648FFF', '#DC267F', '#785EF0', '#FFB000'};
figure
hold on
for iMTF = 1:5
	plot(num_subset, accuracy(iMTF,:), 'Color', colorsMTF{iMTF}, 'linewidth', 2)
end
ylabel('Accuracy')
xlabel('# Neurons in Model')
legend(MTF_names)
grid on 
box off

% Save accuracies
% save(fullfile(base, 'model_comparisons', 'Pop_Rate_Timbre_Subset_MTF.mat'), ...
% 	"accuracy", "std_acc", "CF_groups")
save(fullfile(base, 'model_comparisons', 'Model_Pop_Rate_Timbre_Subset_MTF.mat'), ...
	"accuracy", "std_acc", "CF_groups")
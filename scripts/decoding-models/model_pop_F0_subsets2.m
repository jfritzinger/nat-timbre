%% model_pop_timbre_subsets
clear

%% Load in data
target = 'Bassoon';
%target = 'Oboe';
%target = 'Invariant';

[base, ~, ~, ~] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')
load(fullfile(base, 'model_comparisons', ['Pop_Rate_F0_' target '.mat']), ...
	"pop_rate_F0")
F0s = getF0s(target);

%% Get data into proper matrix

% Get best units
beta_weights = pop_rate_F0.imp.ImportanceMean;
[~,originalpos] = sort(abs(beta_weights), 'descend' );
best_ind = originalpos(1:50);
[~,originalpos] = sort(abs(beta_weights), 'ascend' );
worst_ind = originalpos(1:50);

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
% CF_groups = [0, 14000; 0, 1000; 1000 2000; 2000 3000;3000 4000; ...
% 	4000, 6000;6000 9000;9000 14000];
CF_groups = [0, 14000; 0, 2000; 2000 4000; 4000 8000];
num_subset = [1 5:5:90 100:10:200];

%% Model CF Groups 

iii = 0;
nreps = 1;
accuracy_all = NaN(length(CF_groups), length(num_subset), nreps);
for iCF = 1:length(CF_groups)
	ind = CFs_all > CF_groups(iCF, 1) & CFs_all < CF_groups(iCF, 2);

	beta_weights = pop_rate_F0.imp.ImportanceMean;
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
				T = getF0PopTable(nat_data, target, sesh, F0s, num_data, 'classification');

				% Run model with kfold validation
				[trainedClassifier, accuracy, predictions] = trainClassifierPopRateF0(T, target);
				accuracy_all(iCF, inum, irep) = accuracy;

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
std_acc = std(accuracy_all, [],3);
figure
nexttile
hold on
errorbar(num_subset, accuracy', std_acc'/sqrt(5), 'linewidth', 2)
ylabel('Accuracy')
xlabel('# Neurons in Model')
legend('All', 'CF = 0-2 kHz', 'CF = 2-4 kHz', 'CF = 4-8 kHz')
grid on 
box off

% Save accuracies
save(fullfile(base, 'model_comparisons', ['Pop_Rate_F0_' target '_CF.mat']), ...
	"accuracy", "std_acc", "CF_groups")

%% Model MTF Groups

iii = 0;
accuracy_all = NaN(length(MTF_names), length(num_subset), nreps);
for iMTF = 1:length(MTF_names)

	% Get index 
	ind = ind_all(iMTF,:);

	beta_weights = pop_rate_F0.imp.ImportanceMean;
	[~,best_ind] = sort(abs(beta_weights(ind)), 'descend' );
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
				T = getF0PopTable(nat_data, target, sesh, F0s, num_data);

				% Run model with kfold validation
				[trainedClassifier, accuracy, predictions] = trainClassifierPopRateF0(T, target);
				accuracy_all(iMTF, inum, irep) = accuracy;

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
legend(MTF_names)
grid on 
box off

% Save accuracies
% save(fullfile(base, 'model_comparisons', ['Pop_Rate_F0_' target '_MTF.mat']), ...
% 	"accuracy", "std_acc", "CF_groups")

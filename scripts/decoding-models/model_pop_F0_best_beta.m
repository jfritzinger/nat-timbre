%% model_pop_VS_F0
clear

%% Load in data
%target = 'Bassoon';
target = 'Oboe';

[base, ~, ~, ~] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')
load(fullfile(base, 'model_comparisons', 'Pop_Rate_F0_Oboe3.mat'), ...
	"pop_rate_F0")

%% Get correct output of model

if ismac
	fpath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data';
else
	fpath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\';
end
tuning = readtable(fullfile(fpath, 'Tuning.xlsx')); % Load in tuning

% Get bassoon stimulus
listing = dir(fullfile(fpath, 'waveforms', ['*' target '*.wav']));
files = {listing.name};
note_names = extractBetween(files, 'ff.', '.');
[~, index] = ismember(note_names, tuning.Note);
F0s1 = tuning.Frequency(index);
[F0s, order] = sort(F0s1);

% Find all rows with bassoon and oboe
if strcmp(target, 'Bassoon')
	sesh_all = find(~cellfun(@isempty, {nat_data.bass_rate}));
else
	sesh_all = find(~cellfun(@isempty, {nat_data.oboe_rate}));
end
num_data_all = numel(sesh_all);
CFs_all = [nat_data(sesh_all).CF];
putative = {nat_data(sesh_all).putative};
MTFs = {nat_data(sesh_all).MTF};


% Get indices of interest
beta_weights = pop_rate_F0.imp.ImportanceMean;
[~,originalpos] = sort(abs(beta_weights), 'descend' );
best_ind = originalpos(1:100);
[~,originalpos] = sort(abs(beta_weights), 'ascend' );
worst_ind = originalpos(1:100);

num_neurons = [1:4 5:5:50 60:10:100 1:4 5:5:50 60:10:100];
nmodels = length(num_neurons);
timerVal = tic;
poolobj = parpool();
parfor imodel = 1:nmodels

	for irep = 1:10
		if imodel < nmodels/2+1 % 6 good models
			index = best_ind;
		else % 6 bad models
			index = worst_ind;
		end
		num_index = 1:num_neurons(imodel);
		num_data = num_neurons(imodel);
		sesh = sesh_all(index(num_index));

		% Model including all F0s
		% Find all rows with bassoon in them
		if strcmp(target, 'Bassoon') % Bassoon
			data_mat = NaN(length(F0s)*20, num_data);
			for ii = 1:num_data
				X1 = nat_data(sesh(ii)).bass_raterep';
				X2 = reshape(X1, [], 1);
				data_mat(:,ii) = X2;
			end
		else % Oboe
			data_mat = NaN(length(F0s)*20, num_data);
			for ii = 1:num_data
				X1 = nat_data(sesh(ii)).oboe_raterep';
				X2 = reshape(X1, [], 1);
				data_mat(:,ii) = X2;
			end
		end
		num_data = numel(sesh);

		% Create table for model
		response = reshape(repmat(F0s, 1, 20)', 1, []);
		T = array2table(data_mat);
		T.Response = response';

		if strcmp(target, 'Bassoon')
			[trainedClassifier1, validationAccuracy1, validationPredictions1] = ...
				trainClassifierPopRateF0Bass(T);
		else
			[trainedClassifier1, validationAccuracy1, validationPredictions1] = ...
				trainClassifierPopRateF0Oboe(T);
		end
		C = confusionmat(response, validationPredictions1);
		accuracies(imodel, irep) = sum(diag(C))/sum(C, 'all');

	end
end
timer = toc(timerVal);
fprintf('Models took %0.2g minutes\n', timer/60)
%fprintf('%d/%d, %0.2f%% done!\n', imodel, nmodels, imodel/nmodels*100)
delete(poolobj);

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

save(fullfile(base, 'model_comparisons', 'Pop_Rate_F0_Best_Oboe.mat'), ...
	"accuracies", "num_neurons")


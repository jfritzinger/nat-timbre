%%
clear

%% Load in data and model results

[base, ~, ~, ~] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Pop_Rate_Timbre_All.mat'), ...
	"pop_rate_timbre")
load(fullfile(base, 'model_comparisons', 'Data_NT.mat'), 'nat_data')


%% Set up data matrices and run model

% Get indices of interest
beta_weights = pop_rate_timbre.trainedClassifier.ClassificationSVM.Beta;
[~,originalpos] = sort(abs(beta_weights), 'descend' );
best_ind = originalpos(1:50);
[~,originalpos] = sort(abs(beta_weights), 'ascend' );
worst_ind = originalpos(1:50);

% Get 180 matrix of both bassoon and timbre
has_bass = ~cellfun(@isempty, {nat_data.bass_rate});
has_oboe = ~cellfun(@isempty, {nat_data.oboe_rate});
sesh_all = find(has_bass & has_oboe);

% Find all rows with bassoon and oboe
ind_b = 25:40;
ind_o = [1 3:17];

nmodels = 28;
num_neurons = [1:4 5:5:50 1:4 5:5:50];
for imodel = 1:nmodels

	if imodel < 15 % 6 good models 
		index = best_ind;
	else % 6 bad models 
		index = worst_ind;
	end
	num_index = 1:num_neurons(imodel);
	sesh = sesh_all(index(num_index));
	
	% Model including all F0s
	num_data = numel(sesh);
	CFs = [nat_data(sesh).CF];
	putative = {nat_data(sesh).putative};
	data_mat = NaN(2*20, num_data);
	data_mat2 = NaN(2*20*16, num_data);
	for target = 1:16
		for ii = 1:num_data

			% Arrange data for SVM
			X1 = nat_data(sesh(ii)).bass_raterep(:,ind_b(target));
			X2 = nat_data(sesh(ii)).oboe_raterep(:,ind_o(target));
			X = [X1; X2];
			data_mat(:,ii) = X;
		end
		idx = (1:40) + 40*(target-1);
		data_mat2(idx, :) = data_mat;
	end

	% Put data into table
	T = array2table(data_mat2);
	T.Instrument = repmat([ones(20,1); ones(20, 1)*2], 16, 1);

	[trainedClassifier, accuracy, predictions] = trainClassifierPopRateTimbre(T);
	pop_rate_timbre(imodel).trainedClassifier = trainedClassifier;
	pop_rate_timbre(imodel).accuracy = accuracy;
	pop_rate_timbre(imodel).predictions = predictions;
	pop_rate_timbre(imodel).response = T.Instrument;
	pop_rate_timbre(imodel).T = T;
	pop_rate_timbre(imodel).CFs = CFs;
	pop_rate_timbre(imodel).putative = putative;
	pop_rate_timbre(imodel).sesh = sesh;

	% Print out progress
	fprintf('%d/%d, %0.2f%% done!\n', imodel, nmodels, imodel/nmodels*100)
end

%% Plot accuracy 

figure('Position',[1198,449,293,236])
nneurons = [1:4 5:5:50];
accuracy_bad = [pop_rate_timbre(15:28).accuracy]*100;
plot(nneurons, accuracy_bad);

hold on 
accuracy_good = [pop_rate_timbre(1:14).accuracy]*100;
plot(nneurons, accuracy_good);
xlabel('# Neurons in Model')
ylabel('Accuracy (%)')
grid on

legend('Worst', 'Best')
title('Model ')


%% 



%% Save outputs

% save(fullfile(base, 'model_comparisons', 'Pop_Rate_Timbre_All.mat'), ...
% 	"pop_rate_timbre")

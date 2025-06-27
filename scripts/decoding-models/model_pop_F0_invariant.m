%% model_pop_VS_F0
clear

%% Load in data

[base, ~, ~, ~] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')

%% Get data into proper matrix

% Find all rows with bassoon and oboe
[sesh, num_data] = getF0Sessions(nat_data, 'Invariant');
CFs = [nat_data(sesh).CF];
putative = {nat_data(sesh).putative};

% Model including all F0s 
F0s = getF0s('Invariant');
T = getF0PopTable(nat_data, 'Invariant', sesh, F0s, num_data, 'classification');

%% Run model with kfold validation 

[trainedClassifier, accuracy, predictions] = trainClassifierPopRateF0(T, 'Bassoon');
pop_rate_F0.trainedClassifier = trainedClassifier;
pop_rate_F0.accuracy = accuracy;
pop_rate_F0.predictions = predictions;
pop_rate_F0.response = T.Response;
pop_rate_F0.T = T;
pop_rate_F0.CF = CFs;
pop_rate_F0.putative = putative;
pop_rate_F0.sesh = sesh;
pop_rate_F0.MTF = {nat_data(sesh).MTF};
pop_rate_F0.oboe_rate = [nat_data(sesh).oboe_rate];
pop_rate_F0.oboe_rate_std = [nat_data(sesh).oboe_rate_std];
pop_rate_F0.bass_rate = [nat_data(sesh).bass_rate];
pop_rate_F0.bass_rate_std = [nat_data(sesh).bass_rate_std];

figure
C = confusionmat(T.Response, predictions);
confusionchart(C);
accuracy = sum(diag(C)) / sum(C(:)); % Calculate accuracy
title(['Both ' num2str(accuracy)])
pop_rate_F0.C = C;

% Permutation test for importance 
Mdl = trainedClassifier.ClassificationSVM; % used to be Linear
imp = permutationImportance(Mdl, T, 'Response');
pop_rate_F0.imp = imp;


% Make training and test rows
% ncond = 2;
% ind_test = false(ncond*20, 1); % Preallocate a logical array for 800 elements (40 * 20)
% for istim = 1:ncond
% 	index = randperm(20, 4); % Randomly select 4 indices from 1 to 20
% 	ind_test((istim-1)*20 + index) = true; % Set the selected indices to true
% end
% test_mat = data_mat(ind_test,:);
% train_mat = data_mat(~ind_test,:);
% test_response = [ones(1, 4) repmat(2, 1, 4)]';
% train_response = [ones(1, 16) repmat(2, 1, 16)]';
% T_train = array2table(train_mat);
% T_train.Response = train_response;
% T_test = array2table(test_mat);
% T_test.Response = test_response;

%% Run model 100 times with shuffled data 

ind_b = 25:40;
ind_o = [1 3:17];
rng("shuffle"); % Seed based on current time
nreps = 100;
for imodel = 1:nreps
	data_mat = NaN(2*20, num_data);
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

	% Shuffle data (rows)
	data_mat3 = zeros(640, num_data);
	for ind = 1:640
		row = data_mat2(ind, :);
		shuffled_data = row(randperm(numel(row)));
		data_mat3(ind, :) = shuffled_data;
	end

	% Put data into table
	T_new = array2table(data_mat3);
	T_new.Instrument = repmat([ones(20,1); ones(20, 1)*2], 16, 1);

	% Run model with kfold validation
	[~, shuffled_accuracy(imodel)] = trainClassifierPopRateTimbre(T_new);

	% Print out progress
	fprintf('%d/%d, %0.2f%% done!\n', imodel, nreps, imodel/nreps*100)
end

pop_rate_F0.shuffled_accuracy = shuffled_accuracy;

%% Save model output 

save(fullfile(base, 'model_comparisons', 'Pop_Rate_F0_Invariant.mat'), ...
	"pop_rate_F0")


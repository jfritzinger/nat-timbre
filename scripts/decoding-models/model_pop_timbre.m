%% model_pop_VS_F0
clear

%% Load in data

[base, ~, ~, ~] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Data_NT.mat'), 'nat_data')

%% Get correct output of model

tuning = readtable(fullfile(base, 'Tuning.xlsx')); % Load in tuning
target = 'Bassoon';
listing = dir(fullfile(base, 'waveforms', ['*' target '*.wav']));
files = {listing.name};
note_names = extractBetween(files, 'ff.', '.');
[~, index] = ismember(note_names, tuning.Note);
F0s1 = tuning.Frequency(index);
[F0s, order] = sort(F0s1);

%% Get data into proper matrix
type = 'rate';

% Create array of correct responses
response = [ones(1, 20) repmat(2, 1, 20)]';

% Find all rows with bassoon and oboe
has_bass = ~cellfun(@isempty, {nat_data.bass_rate});
has_oboe = ~cellfun(@isempty, {nat_data.oboe_rate});
sesh = find(has_bass & has_oboe);
num_data = numel(sesh);
ind_b = 25:40;
ind_o = [1 3:17];
CFs = [nat_data(sesh).CF];

% Get putative list
putative = {nat_data(sesh).putative};

% % Sort data 
% for target = 1:16
% 	data_mat = NaN(2*20, num_data);
% 	for ii = 1:num_data
% 
% 		if strcmp(type, 'rate')
% 			X1 = nat_data(sesh(ii)).bass_raterep(:,ind_b(target));
% 			X2 = nat_data(sesh(ii)).oboe_raterep(:,ind_o(target));
% 			X = [X1; X2];
% 		else
% 			X1 = nat_data(sesh(ii)).bass_VSrep(ind_b(target),:)';
% 			X2 = nat_data(sesh(ii)).oboe_VSrep(ind_o(target),:)';
% 			X = [X1; X2];
% 			X(X>=0.99) = NaN;
% 			X(isnan(X)) = 0;
% 		end
% 		data_mat(:,ii) = X;
% 	end
% 
% 	% Create table for model
% 	T = array2table(data_mat);
% 	T.Response = response;
% end

% Model including all F0s 
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

% Put data into table
T = array2table(data_mat2);
T.Instrument = repmat([ones(20,1); ones(20, 1)*2], 16, 1);

%% Run model with kfold validation 

[trainedClassifier, accuracy, predictions] = trainClassifierPopRateTimbre(T);
pop_rate_timbre.trainedClassifier = trainedClassifier;
pop_rate_timbre.accuracy = accuracy;
pop_rate_timbre.predictions = predictions;
pop_rate_timbre.response = T.Instrument;
pop_rate_timbre.T = T;
pop_rate_timbre.CFs = CFs;
pop_rate_timbre.putative = putative;
pop_rate_timbre.sesh = sesh;
pop_rate_timbre.MTF = {nat_data(sesh).MTF};
pop_rate_timbre.oboe_rate = [nat_data(sesh).oboe_rate];
pop_rate_timbre.oboe_rate_std = [nat_data(sesh).oboe_rate_std];
pop_rate_timbre.bass_rate = [nat_data(sesh).bass_rate];
pop_rate_timbre.bass_rate_std = [nat_data(sesh).bass_rate_std];

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


% % This code processes the data into the right shape for testing the
% % model.
% inputTable = T_test;
% predictors = inputTable(:, 1:180);
% response = inputTable.Response;
% classNames = [1; 2];
% 
% % Run model with testing data 
% [yfit,scores] = trainedClassifier.predictFcn(inputTable);

%% Run model 100 times with shuffled data 

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

pop_rate_timbre.shuffled_accuracy = shuffled_accuracy;

%% Save model output 

save(fullfile(base, 'model_comparisons', 'Pop_Rate_Timbre_All.mat'), ...
	"pop_rate_timbre")


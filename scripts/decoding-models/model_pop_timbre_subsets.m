%% model_pop_timbre_subsets

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

% Find all rows with bassoon and oboe
has_bass = ~cellfun(@isempty, {nat_data.bass_rate});
has_oboe = ~cellfun(@isempty, {nat_data.oboe_rate});
sesh = find(has_bass & has_oboe);
num_data = numel(sesh);
ind_b = 25:40;
ind_o = [1 3:17];
CFs = [nat_data(sesh).CF];
putative = {nat_data(sesh).putative};


% EDIT %%%%%%%% 
% Get a subset of the data based on MTF, CF, or something else 

%% Model including all F0s 

% Create array of correct responses
response = [ones(1, 20) repmat(2, 1, 20)]';

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

%% Plot outputs 
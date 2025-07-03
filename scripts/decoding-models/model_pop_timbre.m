%% model_pop_timbre
clear

%% Load in data

[base, ~, ~, ~] = getPathsNT();
%load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')
load(fullfile(base, 'model_comparisons',  'Model_NT.mat'), 'nat_model')
nat_data = nat_model;

%% Get data into proper matrix

% Find all rows with bassoon and oboe
[sesh, num_data] = getTimbreSessions(nat_data);
T = getTimbrePopTable(nat_data, 'Rate', sesh, num_data, 'Model');

%% Run model with kfold validation 

[trainedClassifier, accuracy, predictions] = trainClassifierPopRateTimbre(T);
pop_rate_timbre.trainedClassifier = trainedClassifier;
pop_rate_timbre.accuracy = accuracy;
pop_rate_timbre.predictions = predictions;
pop_rate_timbre.response = T.Instrument;
pop_rate_timbre.T = T;
pop_rate_timbre.CFs = [nat_data(sesh).CF];
pop_rate_timbre.putative = {nat_data(sesh).putative};
pop_rate_timbre.sesh = sesh;
pop_rate_timbre.MTF = {nat_data(sesh).MTF};
pop_rate_timbre.oboe_rate = [nat_data(sesh).oboe_rate];
pop_rate_timbre.oboe_rate_std = [nat_data(sesh).oboe_rate_std];
pop_rate_timbre.bass_rate = [nat_data(sesh).bass_rate];
pop_rate_timbre.bass_rate_std = [nat_data(sesh).bass_rate_std];

% Permutation test for importance 
Mdl = trainedClassifier.ClassificationSVM; % used to be Linear
imp = permutationImportance(Mdl, T, 'Instrument');
pop_rate_timbre.imp = imp;

%% Run model 100 times with shuffled data 

rng("shuffle"); % Seed based on current time
nreps = 100;
shuffled_accuracy = NaN(nreps, 1);
for imodel = 1:nreps

	% Shuffle data (rows)
	data_mat2 = table2array(T(:,1:end-1));
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

% save(fullfile(base, 'model_comparisons', 'Pop_Rate_Timbre_All.mat'), ...
% 	"pop_rate_timbre")
save(fullfile(base, 'model_comparisons', 'Model_Pop_Rate_Timbre_All.mat'), ...
	"pop_rate_timbre")
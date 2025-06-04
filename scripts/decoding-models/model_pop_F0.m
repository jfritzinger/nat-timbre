%% model_pop_VS_F0
clear 

%% Load in data 

[base, ~, ~, ~] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')

%% Get correct output of model 
%target = 'Bassoon';
target = 'Oboe';

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

%% Get data into proper matrix 

% Find all rows with bassoon in them
if strcmp(target, 'Bassoon') % Bassoon
	sesh = find(~cellfun(@isempty, {nat_data.bass_rate}));
	num_data = length(sesh);
	data_mat = NaN(length(F0s)*20, num_data);
	for ii = 1:num_data
		X1 = nat_data(sesh(ii)).bass_raterep';
		X2 = reshape(X1, [], 1);
		data_mat(:,ii) = X2;
	end
else % Oboe 
	sesh = find(~cellfun(@isempty, {nat_data.oboe_rate}));
	num_data = length(sesh);
	data_mat = NaN(length(F0s)*20, num_data);
	for ii = 1:num_data
		X1 = nat_data(sesh(ii)).oboe_raterep';
		X2 = reshape(X1, [], 1);
		data_mat(:,ii) = X2;
	end
end
num_data = numel(sesh);
pop_rate_F0.MTF = {nat_data(sesh).MTF};

% Create array of correct responses
response = reshape(repmat(F0s, 1, 20)', 1, []);

% Create table for model
T = array2table(data_mat);
T.Response = response';

% % Use only first 20 F0s 
% T_new = T(1:400, :);

% % Use only odd F0d
% idx = [];
% for k = 0:19
%     idx = [idx, (1 + 40*k):(20 + 40*k)];
% end
% T_new2 = T(idx,:);


%% Run model 

nrep = 100;
best_accuracy = -inf;
for irep = 1:nrep
	% [trainedClassifier1, validationAccuracy1, validationPredictions1] = ...
	% 	trainClassifierPopRateF0(T, F0s);

	if strcmp(target, 'Bassoon')
		[trainedClassifier1, validationAccuracy1, validationPredictions1] = ...
			trainClassifierPopRateF0Bass(T);
	else
		[trainedClassifier1, validationAccuracy1, validationPredictions1] = ...
			trainClassifierPopRateF0Oboe(T);
	end
	C = confusionmat(response, validationPredictions1);

	accuracy = sum(diag(C))/sum(C, 'all');

	if accuracy > best_accuracy
		trainedClassifier = trainedClassifier1;
		validationAccuracy = validationAccuracy1;
		validationPredictions = validationPredictions1;
		ibest = irep;
		best_accuracy = accuracy;
	end

	% Print out progress
	fprintf('%d/%d, %0.2f%% done!\n', irep, nrep, irep/nrep*100)
end

Mdl = trainedClassifier.ClassificationSVM; % used to be Linear
imp = permutationImportance(Mdl);

%% Save model output 

pop_rate_F0.trainedClassifier = trainedClassifier;
pop_rate_F0.validationPredictions = validationPredictions;
pop_rate_F0.response = response;
pop_rate_F0.T = T;
pop_rate_F0.C = C;
pop_rate_F0.validationAccuracy = validationAccuracy;
pop_rate_F0.putative = {nat_data(sesh).putative};
pop_rate_F0.CF = [nat_data(sesh).CF];
pop_rate_F0.MTF = {nat_data(sesh).MTF};
if strcmp(target, 'Bassoon')
	pop_rate_F0.rate = {nat_data(sesh).bass_rate};
	pop_rate_F0.rate_std = {nat_data(sesh).bass_rate_std};
else
	pop_rate_F0.rate = {nat_data(sesh).oboe_rate};
	pop_rate_F0.rate_std = {nat_data(sesh).oboe_rate_std};
end
pop_rate_F0.imp = imp;

%% Get shuffled accuracy 

rng("shuffle"); % Seed based on current time
nreps = 100;
for imodel = 1:nreps


	% Shuffle data (rows)
	data_mat3 = zeros(length(F0s)*20, num_data);
	for ind = 1:length(F0s)*20
		row = data_mat(ind, :);
		shuffled_data = row(randperm(numel(row)));
		data_mat3(ind, :) = shuffled_data;
	end
	data_mat = data_mat3;

	% Put data into table
	T_new = array2table(data_mat);
	T_new.Response = response';

	% Run model with kfold validation
	if strcmp(target, 'Bassoon')
		[trainedClassifier1, validationAccuracy1, validationPredictions1] = ...
			trainClassifierPopRateF0Bass(T_new);
	else
		[trainedClassifier1, validationAccuracy1, validationPredictions1] = ...
			trainClassifierPopRateF0Oboe(T_new);
	end
	shuffled_accuracy(imodel) = validationAccuracy1;

	% Print out progress
	fprintf('%d/%d, %0.2f%% done!\n', imodel, nreps, imodel/nreps*100)
end
pop_rate_F0.shuffled_accuracy = shuffled_accuracy;

%% Plot outputs 

save(fullfile(base, 'model_comparisons', ...
	['Pop_Rate_F0_' target '3.mat']), "pop_rate_F0")

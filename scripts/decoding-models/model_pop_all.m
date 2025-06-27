%% model_pop_all
clear
timerVal = tic;

%% Load in spreadsheet 

[base, datapath, savepath, ppi] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')

%% Find all rows with bassoon and oboe in them

sesh = [];
for ii = 1:length(nat_data)
	rate = nat_data(ii).bass_rate;
	rate2 = nat_data(ii).oboe_rate;
	if ~isempty(rate) && ~isempty(rate2)
		sesh = [sesh ii];
	end
end
num_data = length(sesh);
ind_b = 25:40;
ind_o = [1 3:17];


%% Set responses 
% 75 * 20 = 1500 responses
tuning = readtable(fullfile(base, 'Tuning.xlsx')); % Load in tuning

% Get bassoon stimulus
target = 'Bassoon';
listing = dir(fullfile(base, 'waveforms', ['*' target '*.wav']));
files = {listing.name};
note_names = extractBetween(files, 'ff.', '.');
[~, index] = ismember(note_names, tuning.Note);
F0s_b = round(tuning.Frequency(index));
[F0s_b, ~] = sort(F0s_b);

% Get bassoon stimulus
target = 'Oboe';
listing = dir(fullfile(base, 'waveforms', ['*' target '*.wav']));
files = {listing.name};
note_names = extractBetween(files, 'ff.', '.');
[~, index] = ismember(note_names, tuning.Note);
F0s_o = round(tuning.Frequency(index));
[F0s_o, ~] = sort(F0s_o);

% Get into 'Response' 
response_b = cell(75,1);
for ii = 1:length(F0s_b)
	response_b{ii} = ['B_' num2str(F0s_b(ii))];
end
for ii = 1:length(F0s_o)
	response_b{ii+length(F0s_b)} = ['O_' num2str(F0s_o(ii))];
end

response = reshape(repmat(response_b, 1, 20)', 1, []);
response = response';

%% Get data 

% Find all rows with bassoon and oboe
has_bass = ~cellfun(@isempty, {nat_data.bass_rate});
has_oboe = ~cellfun(@isempty, {nat_data.oboe_rate});
sesh = find(has_bass & has_oboe);
num_data = numel(sesh);

% Get data 
data_mat1 = NaN(length(F0s_b)*20, num_data);
for ii = 1:num_data
	X1 = nat_data(sesh(ii)).bass_raterep';
	X2 = reshape(X1', [], 1);
	data_mat1(:,ii) = X2;
end

data_mat2 = NaN(length(F0s_o)*20, num_data);
for ii = 1:num_data
	X1 = nat_data(sesh(ii)).oboe_raterep';
	X2 = reshape(X1', [], 1);
	data_mat2(:,ii) = X2;
end

data_mat = [data_mat1; data_mat2];

%% Model! 

T = array2table(data_mat);
T.Response = response;

nrep = 1;
best_accuracy = -inf;
for irep = 1:nrep

	% Model takes about 6 minutes to run
	[trainedClassifier1, validationAccuracy1, validationPredictions1] =...
		trainClassifierAll(T);
	disp(['Model took ' num2str(toc(timerVal)/60) ' minutes'])
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


%% Plot confusion matrix and accuracy 

% Compute classification using all data 
figure
confusionchart(C)

% Calculate accuracy
chart = confusionchart(response,validationPredictions); % Generate confusion chart
confusionMatrix = chart.NormalizedValues; % Get the normalized confusion matrix
accuracy(1) = sum(diag(confusionMatrix)) / sum(confusionMatrix(:)); % Calculate accuracy
title(sprintf('Accuracy = %0.2f%%', accuracy(1)*100))


%% Save output 

pop_rate.trainedClassifier = trainedClassifier;
pop_rate.validationPredictions = validationPredictions;
pop_rate.response = response;
pop_rate.T = T;
pop_rate.C = C;
pop_rate.validationAccuracy = validationAccuracy;
pop_rate.putative = {nat_data(sesh).putative};
pop_rate.CF = [nat_data(sesh).CF];
pop_rate.MTF = {nat_data(sesh).MTF};

save(fullfile(base, 'model_comparisons', 'Pop_Rate2.mat'), "pop_rate")





%% model_pop_all
clear
timerVal = tic;

%% Load in spreadsheet 

[base, datapath, savepath, ppi] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')

%% Set responses 

% Get stimulus
F0s_b = getF0s('Bassoon');
F0s_o = getF0s('Oboe');

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
[sesh, num_data] = getTimbreSessions(nat_data);

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
T = array2table(data_mat);
T.Response = response;

%% Model! 

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





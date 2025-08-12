%% model_pop_all
clear
timerVal = tic;

%% Load in spreadsheet 

[base, datapath, savepath, ppi] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')

%% Model! 

T = getAllPopTable(nat_data, 'Rate');

nrep = 1;
best_accuracy = -inf;
for irep = 1:nrep

	% Separate out training and testing data
	[T_train, T_test] = splitData(T); 

	% Train model
	[trainedClassifier, validationAccuracy, validationPredictions, ~] =...
		trainClassifierAll(T_train);
	disp(['Model took ' num2str(toc(timerVal)/60) ' minutes'])

	% To make predictions with the returned 'trainedClassifier' on new data
	[yfit,scores] = trainedClassifier.predictFcn(T_test);
	C = confusionmat(T_test.Response, yfit);
	accuracy = sum(diag(C))/sum(C, 'all');

	% Print out progress
	fprintf('%d/%d, %0.2f%% done!\n', irep, nrep, irep/nrep*100)
end

% Mdl = trainedClassifier.ClassificationSVM; % used to be Linear
% imp = permutationImportance(Mdl);


%% Plot confusion matrix and accuracy 

% Compute classification using all data 
figure
confusionchart(C)

% Calculate accuracy
chart = confusionchart(T_test.Response,yfit); % Generate confusion chart
confusionMatrix = chart.NormalizedValues; % Get the normalized confusion matrix
title(sprintf('Accuracy = %0.2f%%', accuracy*100))


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

save(fullfile(base, 'model_comparisons', 'Pop_Rate.mat'), "pop_rate")





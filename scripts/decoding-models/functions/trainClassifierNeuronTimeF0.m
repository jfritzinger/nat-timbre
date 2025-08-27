function [trainedClassifier, validationPredictions] = trainClassifierNeuronTimeF0(T, F0s)

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
inputTable = T;
predictorNames = inputTable.Properties.VariableNames(1:end-1);
predictors = inputTable(:, predictorNames);
response = inputTable.Response;
classNames = F0s;

% Train a classifier
% This code specifies all the classifier options and trains the classifier.
template = templateLinear(...
	'Learner', 'SVM', ...
	'Lambda', 'auto', ...
	'BetaTolerance', 0.0001);
classificationLinear = fitcecoc(...
	predictors, ...
	response, ...
	'Learners', template, ...
	'Coding', 'onevsone', ...
	'ClassNames', classNames);

% Create the result struct with predict function
predictorExtractionFcn = @(t) t(:, predictorNames);
classificationLinearPredictFcn = @(x) predict(classificationLinear, x);
trainedClassifier.predictFcn = @(x) classificationLinearPredictFcn(predictorExtractionFcn(x));

% Add additional fields to the result struct
trainedClassifier.RequiredVariables = {'h_all1', 'h_all2', 'h_all3', 'h_all4', 'h_all5', 'h_all6', 'h_all7', 'h_all8', 'h_all9', 'h_all10', 'h_all11', 'h_all12', 'h_all13', 'h_all14', 'h_all15'};
trainedClassifier.ClassificationLinear = classificationLinear;
trainedClassifier.About = 'This struct is a trained model exported from Classification Learner R2025a.';
trainedClassifier.HowToPredict = sprintf('To make predictions on a new table, T, use: \n  [yfit,scores] = c.predictFcn(T) \nreplace ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nThe table, T, must contain the variables returned by: \n  c.RequiredVariables \nVariable formats (e.g. matrix/vector, datatype) must match the original training data. \nAdditional variables are ignored. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appclassification_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
inputTable = T;
predictorNames = inputTable.Properties.VariableNames(1:end-1);
predictors = inputTable(:, predictorNames);
response = inputTable.Response;

% Perform cross-validation
KFolds = 5;
cvp = cvpartition(response, 'KFold', KFolds);

% Initialize the predictions to the proper sizes
validationPredictions = response;
numObservations = size(predictors, 1);
numClasses = length(F0s);
validationScores = NaN(numObservations, numClasses);
for fold = 1:KFolds
	trainingPredictors = predictors(cvp.training(fold), :);
	trainingResponse = response(cvp.training(fold), :);

	% Train a classifier
	% This code specifies all the classifier options and trains the classifier.
	template = templateLinear(...
		'Learner', 'SVM', ...
		'Lambda', 'auto', ...
		'BetaTolerance', 0.0001);
	classificationLinear = fitcecoc(...
		trainingPredictors, ...
		trainingResponse, ...
		'Learners', template, ...
		'Coding', 'onevsone', ...
		'ClassNames', classNames);

	% Create the result struct with predict function
	classificationLinearPredictFcn = @(x) predict(classificationLinear, x);
	validationPredictFcn = @(x) classificationLinearPredictFcn(x);

	% Add additional fields to the result struct

	% Compute validation predictions
	validationPredictors = predictors(cvp.test(fold), :);
	[foldPredictions, foldScores] = validationPredictFcn(validationPredictors);

	% Store predictions in the original order
	validationPredictions(cvp.test(fold), :) = foldPredictions;
	validationScores(cvp.test(fold), :) = foldScores;
end

% Compute validation accuracy
correctPredictions = (validationPredictions == response);
isMissing = isnan(response);
correctPredictions = correctPredictions(~isMissing);
validationAccuracy = sum(correctPredictions)/length(correctPredictions);

end
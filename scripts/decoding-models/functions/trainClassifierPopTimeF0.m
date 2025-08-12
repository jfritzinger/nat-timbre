function [trainedClassifier, validationAccuracy, validationPredictions] ...
	= trainClassifierPopTimeF0(trainingData, F0s)
% [trainedClassifier, validationAccuracy] = trainClassifier(trainingData)
% Returns a trained classifier and its accuracy. This code recreates the
% classification model trained in Classification Learner app. Use the
% generated code to automate training the same model with new data, or to
% learn how to programmatically train models.
%
%  Input:
%      trainingData: A table containing the same predictor and response
%       columns as those imported into the app.
%
%
%  Output:
%      trainedClassifier: A struct containing the trained classifier. The
%       struct contains various fields with information about the trained
%       classifier.
%
%      trainedClassifier.predictFcn: A function to make predictions on new
%       data.
%
%      validationAccuracy: A double representing the validation accuracy as
%       a percentage. In the app, the Models pane displays the validation
%       accuracy for each model.
%
% Use the code to train the model with new data. To retrain your
% classifier, call the function from the command line with your original
% data or new data as the input argument trainingData.
%
% For example, to retrain a classifier trained with the original data set
% T, enter:
%   [trainedClassifier, validationAccuracy] = trainClassifier(T)
%
% To make predictions with the returned 'trainedClassifier' on new data T2,
% use
%   [yfit,scores] = trainedClassifier.predictFcn(T2)
%
% T2 must be a table containing at least the same predictor columns as used
% during training. For details, enter:
%   trainedClassifier.HowToPredict


% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
inputTable = trainingData;
predictors = inputTable(:, 1:end-1);
response = inputTable.response;

% Train a classifier
% This code specifies all the classifier options and trains the classifier.
template = templateLinear(...
	'Learner', 'Logistic', ...
	'Lambda', 'auto', ...
	'BetaTolerance', 0.0001);
classificationLinear = fitcecoc(...
	predictors, ...
	response, ...
	'Learners', template, ...
	'Coding', 'onevsone');

% Create the result struct with predict function
predictorExtractionFcn = @(t) t(:, 1:end-1);
classificationLinearPredictFcn = @(x) predict(classificationLinear, x);
trainedClassifier.predictFcn = @(x) classificationLinearPredictFcn(predictorExtractionFcn(x));

% Add additional fields to the result struct
trainedClassifier.ClassificationLinear = classificationLinear;
trainedClassifier.About = 'This struct is a trained model exported from Classification Learner R2023b.';
trainedClassifier.HowToPredict = sprintf('To make predictions on a new table, T, use: \n  [yfit,scores] = c.predictFcn(T) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nThe table, T, must contain the variables returned by: \n  c.RequiredVariables \nVariable formats (e.g. matrix/vector, datatype) must match the original training data. \nAdditional variables are ignored. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appclassification_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
inputTable = trainingData;
predictors = inputTable(:, 1:end-1);
response = inputTable.response;

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
		'Learner', 'Logistic', ...
		'Lambda', 'auto', ...
		'BetaTolerance', 0.0001);
	classificationLinear = fitcecoc(...
		trainingPredictors, ...
		trainingResponse, ...
		'Learners', template, ...
		'Coding', 'onevsone');

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
correctPredictions = strcmp(validationPredictions, response);
%correctPredictions = (validationPredictions == response);
% = isnan(response);
%correctPredictions = correctPredictions(~isMissing);
validationAccuracy = sum(correctPredictions)/length(correctPredictions);

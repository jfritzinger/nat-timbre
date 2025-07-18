function [trainedClassifier, validationAccuracy, validationPredictions] = ...
	trainClassifierPopRateF0Bass(trainingData, instrument)
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

% Auto-generated by MATLAB on 2025-06-03_042730


% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
inputTable = trainingData;
predictors = inputTable(:, 1:end-1);
response = inputTable.Response;

% Train a classifier
% This code specifies all the classifier options and trains the classifier.
if strcmp(instrument, 'Bassoon')
	template = templateSVM(...
		'KernelFunction', 'gaussian', ...
		'PolynomialOrder', [], ...
		'KernelScale', 993.354569375432, ...
		'BoxConstraint', 0.0474703580635774, ...
		'Standardize', false);
elseif strcmp(instrument, 'Oboe')
	template = templateSVM(...
		'KernelFunction', 'gaussian', ...
		'PolynomialOrder', [], ...
		'KernelScale', 994.1679990352178, ...
		'BoxConstraint', 1.526575074825504, ...
		'Standardize', false);
else
	template = templateSVM(...
		'KernelFunction', 'gaussian', ...
		'PolynomialOrder', [], ...
		'KernelScale', 993.354569375432, ...
		'BoxConstraint', 0.0474703580635774, ...
		'Standardize', false);	
end
classificationSVM = fitcecoc(...
	predictors, ...
	response, ...
	'Learners', template, ...
	'Coding', 'onevsall');

% Create the result struct with predict function
predictorExtractionFcn = @(t) t(:, predictorNames);
svmPredictFcn = @(x) predict(classificationSVM, x);
trainedClassifier.predictFcn = @(x) svmPredictFcn(predictorExtractionFcn(x));

% Add additional fields to the result struct
trainedClassifier.ClassificationSVM = classificationSVM;
trainedClassifier.About = 'This struct is a trained model exported from Classification Learner R2025a.';
trainedClassifier.HowToPredict = sprintf('To make predictions on a new table, T, use: \n  [yfit,scores] = c.predictFcn(T) \nreplace ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nThe table, T, must contain the variables returned by: \n  c.RequiredVariables \nVariable formats (e.g. matrix/vector, datatype) must match the original training data. \nAdditional variables are ignored. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appclassification_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.

% Perform cross-validation
partitionedModel = crossval(trainedClassifier.ClassificationSVM, 'KFold', 5);

% Compute validation predictions
[validationPredictions, validationScores] = kfoldPredict(partitionedModel);

% Compute validation accuracy
validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');

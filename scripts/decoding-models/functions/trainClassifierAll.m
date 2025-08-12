function [trainedClassifier, validationAccuracy, validationPredictions, validationScores] =...
	trainClassifierAll(trainingData)
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
%  Output:
%      trainedClassifier: A struct containing the trained classifier. The
%       struct contains various fields with information about the trained
%       classifier.
%
%      trainedClassifier.predictFcn: A function to make predictions on new
%       data.
%
%      validationAccuracy: A double containing the accuracy as a
%       percentage. In the app, the Models pane displays this overall
%       accuracy score for each model.
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
%   yfit = trainedClassifier.predictFcn(T2)
%
% T2 must be a table containing at least the same predictor columns as used
% during training. For details, enter:
%   trainedClassifier.HowToPredict


% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
inputTable = trainingData;
predictorNames = inputTable.Properties.VariableNames(1:end-1);
predictors = inputTable(:, predictorNames);
response = inputTable.Response;

% Train a classifier
template = templateSVM(...
	'KernelFunction', 'linear', ...
	'PolynomialOrder', [], ...
	'KernelScale', 'auto', ...
	'BoxConstraint', 1, ...
	'Standardize', true);
classificationSVM = fitcecoc(...
	predictors, ...
	response, ...
	'Learners', template, ...
	'Coding', 'onevsone', ...
	'ClassNames', {'B_104'; 'B_110'; 'B_117'; 'B_123'; 'B_131'; 'B_139'; ...
	'B_147'; 'B_156'; 'B_165'; 'B_175'; 'B_185'; 'B_196'; 'B_208'; 'B_220'; ...
	'B_233'; 'B_262'; 'B_277'; 'B_294'; 'B_311'; 'B_330'; 'B_349'; 'B_370'; ...
	'B_392'; 'B_415'; 'B_440'; 'B_466'; 'B_494'; 'B_523'; 'B_554'; 'B_58'; ...
	'B_587'; 'B_62'; 'B_65'; 'B_69'; 'B_73'; 'B_78'; 'B_82'; 'B_87'; 'B_93'; ...
	'B_98'; 'O_1047'; 'O_1109'; 'O_1175'; 'O_1245'; 'O_1319'; 'O_1397'; ...
	'O_1480'; 'O_1568'; 'O_1661'; 'O_233'; 'O_247'; 'O_262'; 'O_277'; ...
	'O_294'; 'O_311'; 'O_330'; 'O_349'; 'O_370'; 'O_392'; 'O_415'; ...
	'O_440'; 'O_466'; 'O_494'; 'O_523'; 'O_554'; 'O_587'; 'O_622'; ...
	'O_659'; 'O_698'; 'O_740'; 'O_784'; 'O_831'; 'O_880'; 'O_932'; 'O_988'});


% Create the result struct with predict function
predictorExtractionFcn = @(t) t(:, predictorNames);
svmPredictFcn = @(x) predict(classificationSVM, x);
trainedClassifier.predictFcn = @(x) svmPredictFcn(predictorExtractionFcn(x));

% Add additional fields to the result struct
trainedClassifier.ClassificationSVM = classificationSVM;
trainedClassifier.HowToPredict = sprintf(['To make predictions on a new table,' ...
	'     T, use: \n  yfit = c.predictFcn(T) \nreplacing ''c'' with the ' ...
	'    name of the variable that is this struct, e.g. ''trainedModel''.' ...
	'     \n \nThe table, T, must contain the variables returned by: \n  ' ...
	'    c.RequiredVariables \nVariable formats (e.g. matrix/vector, datatype) ' ...
	'    must match the original training data. \nAdditional variables are ' ...
	'    ignored. \n \nFor more information, see ' ...
	'    <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map'')' ...
	'    , ''appclassification_exportmodeltoworkspace'')">' ...
	'    How to predict using an exported model</a>.']);

% Perform cross-validation
partitionedModel = crossval(trainedClassifier.ClassificationSVM, 'KFold', 5);

% Compute validation predictions
[validationPredictions, validationScores] = kfoldPredict(partitionedModel);

% Compute validation accuracy
validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');

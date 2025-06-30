function [validationPredictions] = trainClassifierNeuronTimeF0(T, F0s)


response = T.Response;
%predictors = T(1:end-1,:);
predictors = T(:,1:end-1);

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

	% Compute validation predictions
	validationPredictors = predictors(cvp.test(fold), :);
	[foldPredictions, foldScores] = validationPredictFcn(validationPredictors);

	% Store predictions in the original order
	validationPredictions(cvp.test(fold), :) = foldPredictions;
	validationScores(cvp.test(fold), :) = foldScores;
end

end
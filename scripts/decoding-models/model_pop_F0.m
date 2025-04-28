%% model_pop_VS_F0
clear 

%% Load in data 

filepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/model_comparisons';
load(fullfile(filepath, 'Data_NT.mat'), 'nat_data')

%% Get correct output of model 
target = 'Bassoon';
%target = 'Oboe';

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
%F0s = log10(F0s);

%% Get data into proper matrix 
type = 'rate'; % VS

% Find all rows with bassoon in them
if strcmp(target, 'Bassoon')
	sesh = find(~cellfun(@isempty, {nat_data.bass_rate}));
else
	sesh = find(~cellfun(@isempty, {nat_data.oboe_rate}));
end
num_data = numel(sesh);

data_mat = NaN(length(F0s)*20, num_data);
for ii = 1:num_data

	if strcmp(type, 'rate')
		if strcmp(target, 'Bassoon')
			X1 = nat_data(sesh(ii)).bass_raterep';
		else
			X1 = nat_data(sesh(ii)).oboe_raterep';
		end
	else
		if strcmp(target, 'Bassoon')
			X1 = nat_data(sesh(ii)).bass_VSrep;
		else
			X1 = nat_data(sesh(ii)).oboe_VSrep;
		end
	end

	X2 = reshape(X1, [], 1);
	data_mat(:,ii) = X2;
end

% Create array of correct responses
response = reshape(repmat(F0s, 1, 20)', 1, []);

% Create table for model
T = array2table(data_mat);
T.Response = response';

%% Run model 

nrep = 50;
for irep = 1:nrep
	% Extract predictors and response
	% This code processes the data into the right shape for training the
	% model.
	inputTable = T;
	predictorNames = {'data_mat1', 'data_mat2', 'data_mat3', 'data_mat4', 'data_mat5', 'data_mat6', 'data_mat7', 'data_mat8', 'data_mat9', 'data_mat10', 'data_mat11', 'data_mat12', 'data_mat13', 'data_mat14', 'data_mat15', 'data_mat16', 'data_mat17', 'data_mat18', 'data_mat19', 'data_mat20', 'data_mat21', 'data_mat22', 'data_mat23', 'data_mat24', 'data_mat25', 'data_mat26', 'data_mat27', 'data_mat28', 'data_mat29', 'data_mat30', 'data_mat31', 'data_mat32', 'data_mat33', 'data_mat34', 'data_mat35', 'data_mat36', 'data_mat37', 'data_mat38', 'data_mat39', 'data_mat40', 'data_mat41', 'data_mat42', 'data_mat43', 'data_mat44', 'data_mat45', 'data_mat46', 'data_mat47', 'data_mat48', 'data_mat49', 'data_mat50', 'data_mat51', 'data_mat52', 'data_mat53', 'data_mat54', 'data_mat55', 'data_mat56', 'data_mat57', 'data_mat58', 'data_mat59', 'data_mat60', 'data_mat61', 'data_mat62', 'data_mat63', 'data_mat64', 'data_mat65', 'data_mat66', 'data_mat67', 'data_mat68', 'data_mat69', 'data_mat70', 'data_mat71', 'data_mat72', 'data_mat73', 'data_mat74', 'data_mat75', 'data_mat76', 'data_mat77', 'data_mat78', 'data_mat79', 'data_mat80', 'data_mat81', 'data_mat82', 'data_mat83', 'data_mat84', 'data_mat85', 'data_mat86', 'data_mat87', 'data_mat88', 'data_mat89', 'data_mat90', 'data_mat91', 'data_mat92', 'data_mat93', 'data_mat94', 'data_mat95', 'data_mat96', 'data_mat97', 'data_mat98', 'data_mat99', 'data_mat100', 'data_mat101', 'data_mat102', 'data_mat103', 'data_mat104', 'data_mat105', 'data_mat106', 'data_mat107', 'data_mat108', 'data_mat109', 'data_mat110', 'data_mat111', 'data_mat112', 'data_mat113', 'data_mat114', 'data_mat115', 'data_mat116', 'data_mat117', 'data_mat118', 'data_mat119', 'data_mat120', 'data_mat121', 'data_mat122', 'data_mat123', 'data_mat124', 'data_mat125', 'data_mat126', 'data_mat127', 'data_mat128', 'data_mat129', 'data_mat130', 'data_mat131', 'data_mat132', 'data_mat133', 'data_mat134', 'data_mat135', 'data_mat136', 'data_mat137', 'data_mat138', 'data_mat139', 'data_mat140', 'data_mat141', 'data_mat142', 'data_mat143', 'data_mat144', 'data_mat145', 'data_mat146', 'data_mat147', 'data_mat148', 'data_mat149', 'data_mat150', 'data_mat151', 'data_mat152', 'data_mat153', 'data_mat154', 'data_mat155', 'data_mat156', 'data_mat157', 'data_mat158', 'data_mat159', 'data_mat160', 'data_mat161', 'data_mat162', 'data_mat163', 'data_mat164', 'data_mat165', 'data_mat166', 'data_mat167', 'data_mat168', 'data_mat169', 'data_mat170', 'data_mat171', 'data_mat172', 'data_mat173', 'data_mat174', 'data_mat175', 'data_mat176', 'data_mat177', 'data_mat178', 'data_mat179', 'data_mat180', 'data_mat181', 'data_mat182', 'data_mat183', 'data_mat184', 'data_mat185', 'data_mat186', 'data_mat187', 'data_mat188', 'data_mat189', 'data_mat190', 'data_mat191', 'data_mat192', 'data_mat193', 'data_mat194', 'data_mat195', 'data_mat196', 'data_mat197', 'data_mat198', 'data_mat199', 'data_mat200', 'data_mat201', 'data_mat202', 'data_mat203', 'data_mat204', 'data_mat205', 'data_mat206', 'data_mat207', 'data_mat208', 'data_mat209', 'data_mat210', 'data_mat211', 'data_mat212', 'data_mat213'};
	predictors = inputTable(:, predictorNames);
	response = inputTable.Response;
	isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false];
	classNames = [58.27; 61.74; 65.41; 69.3; 73.42; 77.78; 82.41; 87.31; 92.5; 98; 103.83; 110; 116.54; 123.47; 130.81; 138.59; 146.83; 155.56; 164.81; 174.61; 185; 196; 207.65; 220; 233.08; 261.63; 277.18; 293.66; 311.13; 329.63; 349.23; 369.99; 392; 415.3; 440; 466.16; 493.88; 523.25; 554.37; 587.33];

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
	trainedClassifier.RequiredVariables = {'data_mat1', 'data_mat2', 'data_mat3', 'data_mat4', 'data_mat5', 'data_mat6', 'data_mat7', 'data_mat8', 'data_mat9', 'data_mat10', 'data_mat11', 'data_mat12', 'data_mat13', 'data_mat14', 'data_mat15', 'data_mat16', 'data_mat17', 'data_mat18', 'data_mat19', 'data_mat20', 'data_mat21', 'data_mat22', 'data_mat23', 'data_mat24', 'data_mat25', 'data_mat26', 'data_mat27', 'data_mat28', 'data_mat29', 'data_mat30', 'data_mat31', 'data_mat32', 'data_mat33', 'data_mat34', 'data_mat35', 'data_mat36', 'data_mat37', 'data_mat38', 'data_mat39', 'data_mat40', 'data_mat41', 'data_mat42', 'data_mat43', 'data_mat44', 'data_mat45', 'data_mat46', 'data_mat47', 'data_mat48', 'data_mat49', 'data_mat50', 'data_mat51', 'data_mat52', 'data_mat53', 'data_mat54', 'data_mat55', 'data_mat56', 'data_mat57', 'data_mat58', 'data_mat59', 'data_mat60', 'data_mat61', 'data_mat62', 'data_mat63', 'data_mat64', 'data_mat65', 'data_mat66', 'data_mat67', 'data_mat68', 'data_mat69', 'data_mat70', 'data_mat71', 'data_mat72', 'data_mat73', 'data_mat74', 'data_mat75', 'data_mat76', 'data_mat77', 'data_mat78', 'data_mat79', 'data_mat80', 'data_mat81', 'data_mat82', 'data_mat83', 'data_mat84', 'data_mat85', 'data_mat86', 'data_mat87', 'data_mat88', 'data_mat89', 'data_mat90', 'data_mat91', 'data_mat92', 'data_mat93', 'data_mat94', 'data_mat95', 'data_mat96', 'data_mat97', 'data_mat98', 'data_mat99', 'data_mat100', 'data_mat101', 'data_mat102', 'data_mat103', 'data_mat104', 'data_mat105', 'data_mat106', 'data_mat107', 'data_mat108', 'data_mat109', 'data_mat110', 'data_mat111', 'data_mat112', 'data_mat113', 'data_mat114', 'data_mat115', 'data_mat116', 'data_mat117', 'data_mat118', 'data_mat119', 'data_mat120', 'data_mat121', 'data_mat122', 'data_mat123', 'data_mat124', 'data_mat125', 'data_mat126', 'data_mat127', 'data_mat128', 'data_mat129', 'data_mat130', 'data_mat131', 'data_mat132', 'data_mat133', 'data_mat134', 'data_mat135', 'data_mat136', 'data_mat137', 'data_mat138', 'data_mat139', 'data_mat140', 'data_mat141', 'data_mat142', 'data_mat143', 'data_mat144', 'data_mat145', 'data_mat146', 'data_mat147', 'data_mat148', 'data_mat149', 'data_mat150', 'data_mat151', 'data_mat152', 'data_mat153', 'data_mat154', 'data_mat155', 'data_mat156', 'data_mat157', 'data_mat158', 'data_mat159', 'data_mat160', 'data_mat161', 'data_mat162', 'data_mat163', 'data_mat164', 'data_mat165', 'data_mat166', 'data_mat167', 'data_mat168', 'data_mat169', 'data_mat170', 'data_mat171', 'data_mat172', 'data_mat173', 'data_mat174', 'data_mat175', 'data_mat176', 'data_mat177', 'data_mat178', 'data_mat179', 'data_mat180', 'data_mat181', 'data_mat182', 'data_mat183', 'data_mat184', 'data_mat185', 'data_mat186', 'data_mat187', 'data_mat188', 'data_mat189', 'data_mat190', 'data_mat191', 'data_mat192', 'data_mat193', 'data_mat194', 'data_mat195', 'data_mat196', 'data_mat197', 'data_mat198', 'data_mat199', 'data_mat200', 'data_mat201', 'data_mat202', 'data_mat203', 'data_mat204', 'data_mat205', 'data_mat206', 'data_mat207', 'data_mat208', 'data_mat209', 'data_mat210', 'data_mat211', 'data_mat212', 'data_mat213'};
	trainedClassifier.ClassificationLinear = classificationLinear;
	trainedClassifier.About = 'This struct is a trained model exported from Classification Learner R2023b.';
	trainedClassifier.HowToPredict = sprintf('To make predictions on a new table, T, use: \n  [yfit,scores] = c.predictFcn(T) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nThe table, T, must contain the variables returned by: \n  c.RequiredVariables \nVariable formats (e.g. matrix/vector, datatype) must match the original training data. \nAdditional variables are ignored. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appclassification_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

	% Extract predictors and response
	% This code processes the data into the right shape for training the
	% model.
	inputTable = T;
	predictorNames = {'data_mat1', 'data_mat2', 'data_mat3', 'data_mat4', 'data_mat5', 'data_mat6', 'data_mat7', 'data_mat8', 'data_mat9', 'data_mat10', 'data_mat11', 'data_mat12', 'data_mat13', 'data_mat14', 'data_mat15', 'data_mat16', 'data_mat17', 'data_mat18', 'data_mat19', 'data_mat20', 'data_mat21', 'data_mat22', 'data_mat23', 'data_mat24', 'data_mat25', 'data_mat26', 'data_mat27', 'data_mat28', 'data_mat29', 'data_mat30', 'data_mat31', 'data_mat32', 'data_mat33', 'data_mat34', 'data_mat35', 'data_mat36', 'data_mat37', 'data_mat38', 'data_mat39', 'data_mat40', 'data_mat41', 'data_mat42', 'data_mat43', 'data_mat44', 'data_mat45', 'data_mat46', 'data_mat47', 'data_mat48', 'data_mat49', 'data_mat50', 'data_mat51', 'data_mat52', 'data_mat53', 'data_mat54', 'data_mat55', 'data_mat56', 'data_mat57', 'data_mat58', 'data_mat59', 'data_mat60', 'data_mat61', 'data_mat62', 'data_mat63', 'data_mat64', 'data_mat65', 'data_mat66', 'data_mat67', 'data_mat68', 'data_mat69', 'data_mat70', 'data_mat71', 'data_mat72', 'data_mat73', 'data_mat74', 'data_mat75', 'data_mat76', 'data_mat77', 'data_mat78', 'data_mat79', 'data_mat80', 'data_mat81', 'data_mat82', 'data_mat83', 'data_mat84', 'data_mat85', 'data_mat86', 'data_mat87', 'data_mat88', 'data_mat89', 'data_mat90', 'data_mat91', 'data_mat92', 'data_mat93', 'data_mat94', 'data_mat95', 'data_mat96', 'data_mat97', 'data_mat98', 'data_mat99', 'data_mat100', 'data_mat101', 'data_mat102', 'data_mat103', 'data_mat104', 'data_mat105', 'data_mat106', 'data_mat107', 'data_mat108', 'data_mat109', 'data_mat110', 'data_mat111', 'data_mat112', 'data_mat113', 'data_mat114', 'data_mat115', 'data_mat116', 'data_mat117', 'data_mat118', 'data_mat119', 'data_mat120', 'data_mat121', 'data_mat122', 'data_mat123', 'data_mat124', 'data_mat125', 'data_mat126', 'data_mat127', 'data_mat128', 'data_mat129', 'data_mat130', 'data_mat131', 'data_mat132', 'data_mat133', 'data_mat134', 'data_mat135', 'data_mat136', 'data_mat137', 'data_mat138', 'data_mat139', 'data_mat140', 'data_mat141', 'data_mat142', 'data_mat143', 'data_mat144', 'data_mat145', 'data_mat146', 'data_mat147', 'data_mat148', 'data_mat149', 'data_mat150', 'data_mat151', 'data_mat152', 'data_mat153', 'data_mat154', 'data_mat155', 'data_mat156', 'data_mat157', 'data_mat158', 'data_mat159', 'data_mat160', 'data_mat161', 'data_mat162', 'data_mat163', 'data_mat164', 'data_mat165', 'data_mat166', 'data_mat167', 'data_mat168', 'data_mat169', 'data_mat170', 'data_mat171', 'data_mat172', 'data_mat173', 'data_mat174', 'data_mat175', 'data_mat176', 'data_mat177', 'data_mat178', 'data_mat179', 'data_mat180', 'data_mat181', 'data_mat182', 'data_mat183', 'data_mat184', 'data_mat185', 'data_mat186', 'data_mat187', 'data_mat188', 'data_mat189', 'data_mat190', 'data_mat191', 'data_mat192', 'data_mat193', 'data_mat194', 'data_mat195', 'data_mat196', 'data_mat197', 'data_mat198', 'data_mat199', 'data_mat200', 'data_mat201', 'data_mat202', 'data_mat203', 'data_mat204', 'data_mat205', 'data_mat206', 'data_mat207', 'data_mat208', 'data_mat209', 'data_mat210', 'data_mat211', 'data_mat212', 'data_mat213'};
	predictors = inputTable(:, predictorNames);
	response = inputTable.Response;
	isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false];
	classNames = [58.27; 61.74; 65.41; 69.3; 73.42; 77.78; 82.41; 87.31; 92.5; 98; 103.83; 110; 116.54; 123.47; 130.81; 138.59; 146.83; 155.56; 164.81; 174.61; 185; 196; 207.65; 220; 233.08; 261.63; 277.18; 293.66; 311.13; 329.63; 349.23; 369.99; 392; 415.3; 440; 466.16; 493.88; 523.25; 554.37; 587.33];

	% Perform cross-validation
	KFolds = 5;
	cvp = cvpartition(response, 'KFold', KFolds);
	% Initialize the predictions to the proper sizes
	validationPredictions = response;
	numObservations = size(predictors, 1);
	numClasses = 40;
	validationScores = NaN(numObservations, numClasses);
	for fold = 1:KFolds
		trainingPredictors = predictors(cvp.training(fold), :);
		trainingResponse = response(cvp.training(fold), :);
		foldIsCategoricalPredictor = isCategoricalPredictor;

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
	validationAccuracy(irep) = sum(correctPredictions)/length(correctPredictions);

	fprintf('Done: %d/%d\n', irep, nrep)
end

%% Plot many model reps 

figure
swarmchart(ones(50,1), validationAccuracy, 'filled')
hold on
boxplot(validationAccuracy)
ylabel('Model Accuracy')
title('50 reps of model')


%% Plot confusion matrix 

% Compute classification using all data 
figure
C = confusionmat(response, validationPredictions);
confusionchart(C)

% Calculate accuracy
chart = confusionchart(response,validationPredictions); % Generate confusion chart
confusionMatrix = chart.NormalizedValues; % Get the normalized confusion matrix
accuracy(1) = sum(diag(confusionMatrix)) / sum(confusionMatrix(:)); % Calculate accuracy
title(sprintf('Accuracy = %0.2f%%', accuracy(1)*100))

% Accurate % 202 / 800
num_acc =  sum(diag(confusionMatrix));

% Accuracy within +/- 1 semitone % 13 / 800
num_acc1 =  sum(diag(confusionMatrix, 1)) + sum(diag(confusionMatrix, -1));

% Accuracy within +/- 2 semitones % 319 / 800
num_acc2 =  sum(diag(confusionMatrix, 2)) + sum(diag(confusionMatrix, -2));

% Accuracy within +/- 3 semitones % 319 / 800
num_acc3 =  sum(diag(confusionMatrix, 3)) + sum(diag(confusionMatrix, -3));

% Accuracy within +/- 4 semitones % 319 / 800
num_acc4 =  sum(diag(confusionMatrix, 4)) + sum(diag(confusionMatrix, -4));

%% Plot model 

figure

nexttile
beta_weights = classificationLinearPredictFcn.Beta;
plot(1:180, beta_weights)
hold on 
yline(0)
xlabel('Neuron #')
ylabel('Beta Weight')

nexttile
scatter(CFs/1000, beta_weights)
hold on 
yline(0)
xticks([0.1 0.2 0.5 1 2 5 10])
xlim([0.2 10])
xlabel('CF')
ylabel('Beta Weight')
set(gca, 'xscale', 'log')

nexttile
scatter(CFs/1000, abs(beta_weights))
hold on 
yline(0)
xticks([0.1 0.2 0.5 1 2 5 10])
xlim([0.2 10])
xlabel('CF')
ylabel('abs(Beta Weight)')
set(gca, 'xscale', 'log')


%% Run model
% 
% nrep = 5;
% for irep = 1:nrep
% 
% 	% Take out test data
% 	ind_test = false(length(F0s)*20, 1); % Preallocate a logical array for 800 elements (40 * 20)
% 	for istim = 1:length(F0s)
% 		index = randperm(20, 4); % Randomly select 4 indices from 1 to 20
% 		ind_test((istim-1)*20 + index) = true; % Set the selected indices to true
% 	end
% 
% 	% Make training and test rows
% 	test_mat = data_mat(ind_test);
% 	train_mat = data_mat(~ind_test);
% 
% 	% Train model on training data
% 	% mdl = fitrlinear(train_mat, response,'Regularization','lasso',...
% 	% 	'Solver','sparsa', 'Lambda','auto');
% 	%mdl = fitlm(train_mat, response);
% 	% mdl = fitrlinear(train_mat, response, 'Regularization', 'lasso', ...
% 	% 	'Lambda', 0.01, 'KFold', 5); % Linear discriminant analysis
% 	mdl = fitcecoc(train_mat, response, ...
% 		'OptimizeHyperparameters','auto'); % Multiclass SVM using Error-Correcting Output Codes
% 
% 	% If the loss is very high (close to random guessing), 
% 	% revisit your data preprocessing steps.	
% 	%cvMdl = fitcecoc(train_mat, response, 'KFold', 5);
% 	% cvLoss = kfoldLoss(mdl); % Cross-validation loss
% 	% disp(['Cross-validation loss: ', num2str(cvLoss)]);
% 
% 	% Evaluate
% 	% trainedModel = mdl.Trained{1}; % Access the first fold's trained model
% 	output = predict(mdl, test_mat);
% 	output_avg = mean(reshape(output, length(F0s), 4), 2);
% 
% 	% Subtract real - predicted to get accuracy
% 	accuracy2(irep,:) = sum(output == response_test') / length(response_test);
% 	accuracy(irep,:) = F0s - output_avg;
% 	alloutput_avg(irep, :) = output_avg;
% 
% 	% Correlation coefficient
% 	R_temp = corrcoef(response_test, output);
% 	R(irep) = R_temp(1, 2);
% 	R2(irep) = R_temp(1,2)^2;
% 
% end
% 
% %% Plot outputs 
% 
% % Confusion matrix for best model 
% [best_R2, best_ind] = max(R2);
% 
% output_hz = 10.^output_avg;
% figure('Position',[231,839,1016,351])
% nexttile
% scatter(F0s1, output_hz)
% hold on
% plot([1 2000], [1 2000], 'k')
% xlim([min(F0s1), max(F0s1)])
% ylim([min(F0s1), max(F0s1)])
% set(gca, 'XScale', 'log', 'YScale', 'log')
% xlabel('Actual F0 (Hz)')
% ylabel('Predicted F0 (Hz)')
% 
% % Plot accuracy 
% nexttile
% scatter(F0s, accuracy, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.5)
% %set(gca, 'XScale', 'log')
% 
% 
% nexttile
% histogram(accuracy2*100)
% xlim([0 50])
%% model_pop_VS_F0
clear

%% Load in data

filepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/model_comparisons';
load(fullfile(filepath, 'Data_NT.mat'), 'nat_data')

%% Get correct output of model

if ismac
	fpath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data';
else
	fpath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\';
end
tuning = readtable(fullfile(fpath, 'Tuning.xlsx')); % Load in tuning

% Get bassoon stimulus
target = 'Bassoon';
listing = dir(fullfile(fpath, 'waveforms', ['*' target '*.wav']));
files = {listing.name};
note_names = extractBetween(files, 'ff.', '.');
[~, index] = ismember(note_names, tuning.Note);
F0s1 = tuning.Frequency(index);
[F0s, order] = sort(F0s1);

%% Get data into proper matrix
type = 'rate'; % VS

% Create array of correct responses
response = [ones(1, 20) repmat(2, 1, 20)]';

% Find all rows with bassoon and oboe
has_bass = ~cellfun(@isempty, {nat_data.bass_rate});
has_oboe = ~cellfun(@isempty, {nat_data.oboe_rate});
sesh = find(has_bass & has_oboe);
num_data = numel(sesh);
ind_b = 25:40;
ind_o = [1 3:17];
CFs = [nat_data(sesh).CF];

% Model including all F0s 
data_mat = NaN(2*20, num_data);
for target = 1:16
	for ii = 1:num_data

		% Arrange data for SVM
		X1 = nat_data(sesh(ii)).bass_raterep(:,ind_b(target));
		X2 = nat_data(sesh(ii)).oboe_raterep(:,ind_o(target));
		X = [X1; X2];
		data_mat(:,ii) = X;

	end
	idx = (1:40) + 40*(target-1);
	data_mat2(idx, :) = data_mat;

end

% Put data into table
T_new = array2table(data_mat2);
T_new.Instrument = repmat([ones(20,1); ones(20, 1)*2], 16, 1);

%% Model 

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
inputTable = T_new;
predictorNames = {'data_mat21', 'data_mat22', 'data_mat23', 'data_mat24', 'data_mat25', 'data_mat26', 'data_mat27', 'data_mat28', 'data_mat29', 'data_mat210', 'data_mat211', 'data_mat212', 'data_mat213', 'data_mat214', 'data_mat215', 'data_mat216', 'data_mat217', 'data_mat218', 'data_mat219', 'data_mat220', 'data_mat221', 'data_mat222', 'data_mat223', 'data_mat224', 'data_mat225', 'data_mat226', 'data_mat227', 'data_mat228', 'data_mat229', 'data_mat230', 'data_mat231', 'data_mat232', 'data_mat233', 'data_mat234', 'data_mat235', 'data_mat236', 'data_mat237', 'data_mat238', 'data_mat239', 'data_mat240', 'data_mat241', 'data_mat242', 'data_mat243', 'data_mat244', 'data_mat245', 'data_mat246', 'data_mat247', 'data_mat248', 'data_mat249', 'data_mat250', 'data_mat251', 'data_mat252', 'data_mat253', 'data_mat254', 'data_mat255', 'data_mat256', 'data_mat257', 'data_mat258', 'data_mat259', 'data_mat260', 'data_mat261', 'data_mat262', 'data_mat263', 'data_mat264', 'data_mat265', 'data_mat266', 'data_mat267', 'data_mat268', 'data_mat269', 'data_mat270', 'data_mat271', 'data_mat272', 'data_mat273', 'data_mat274', 'data_mat275', 'data_mat276', 'data_mat277', 'data_mat278', 'data_mat279', 'data_mat280', 'data_mat281', 'data_mat282', 'data_mat283', 'data_mat284', 'data_mat285', 'data_mat286', 'data_mat287', 'data_mat288', 'data_mat289', 'data_mat290', 'data_mat291', 'data_mat292', 'data_mat293', 'data_mat294', 'data_mat295', 'data_mat296', 'data_mat297', 'data_mat298', 'data_mat299', 'data_mat2100', 'data_mat2101', 'data_mat2102', 'data_mat2103', 'data_mat2104', 'data_mat2105', 'data_mat2106', 'data_mat2107', 'data_mat2108', 'data_mat2109', 'data_mat2110', 'data_mat2111', 'data_mat2112', 'data_mat2113', 'data_mat2114', 'data_mat2115', 'data_mat2116', 'data_mat2117', 'data_mat2118', 'data_mat2119', 'data_mat2120', 'data_mat2121', 'data_mat2122', 'data_mat2123', 'data_mat2124', 'data_mat2125', 'data_mat2126', 'data_mat2127', 'data_mat2128', 'data_mat2129', 'data_mat2130', 'data_mat2131', 'data_mat2132', 'data_mat2133', 'data_mat2134', 'data_mat2135', 'data_mat2136', 'data_mat2137', 'data_mat2138', 'data_mat2139', 'data_mat2140', 'data_mat2141', 'data_mat2142', 'data_mat2143', 'data_mat2144', 'data_mat2145', 'data_mat2146', 'data_mat2147', 'data_mat2148', 'data_mat2149', 'data_mat2150', 'data_mat2151', 'data_mat2152', 'data_mat2153', 'data_mat2154', 'data_mat2155', 'data_mat2156', 'data_mat2157', 'data_mat2158', 'data_mat2159', 'data_mat2160', 'data_mat2161', 'data_mat2162', 'data_mat2163', 'data_mat2164', 'data_mat2165', 'data_mat2166', 'data_mat2167', 'data_mat2168', 'data_mat2169', 'data_mat2170', 'data_mat2171', 'data_mat2172', 'data_mat2173', 'data_mat2174', 'data_mat2175', 'data_mat2176', 'data_mat2177', 'data_mat2178', 'data_mat2179', 'data_mat2180'};
predictors = inputTable(:, predictorNames);
response = inputTable.Instrument;
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false];
classNames = [1; 2];

% Train a classifier
% This code specifies all the classifier options and trains the classifier.
classificationLinear = fitclinear(...
	predictors, ...
	response, ...
	'Learner', 'SVM', ...
	'Lambda', 'auto', ...
	'BetaTolerance', 0.0001, ...
	'ClassNames', classNames);

% Create the result struct with predict function
predictorExtractionFcn = @(t) t(:, predictorNames);
classificationLinearPredictFcn = @(x) predict(classificationLinear, x);
trainedClassifier.predictFcn = @(x) classificationLinearPredictFcn(predictorExtractionFcn(x));

% Add additional fields to the result struct
trainedClassifier.RequiredVariables = {'data_mat21', 'data_mat22', 'data_mat23', 'data_mat24', 'data_mat25', 'data_mat26', 'data_mat27', 'data_mat28', 'data_mat29', 'data_mat210', 'data_mat211', 'data_mat212', 'data_mat213', 'data_mat214', 'data_mat215', 'data_mat216', 'data_mat217', 'data_mat218', 'data_mat219', 'data_mat220', 'data_mat221', 'data_mat222', 'data_mat223', 'data_mat224', 'data_mat225', 'data_mat226', 'data_mat227', 'data_mat228', 'data_mat229', 'data_mat230', 'data_mat231', 'data_mat232', 'data_mat233', 'data_mat234', 'data_mat235', 'data_mat236', 'data_mat237', 'data_mat238', 'data_mat239', 'data_mat240', 'data_mat241', 'data_mat242', 'data_mat243', 'data_mat244', 'data_mat245', 'data_mat246', 'data_mat247', 'data_mat248', 'data_mat249', 'data_mat250', 'data_mat251', 'data_mat252', 'data_mat253', 'data_mat254', 'data_mat255', 'data_mat256', 'data_mat257', 'data_mat258', 'data_mat259', 'data_mat260', 'data_mat261', 'data_mat262', 'data_mat263', 'data_mat264', 'data_mat265', 'data_mat266', 'data_mat267', 'data_mat268', 'data_mat269', 'data_mat270', 'data_mat271', 'data_mat272', 'data_mat273', 'data_mat274', 'data_mat275', 'data_mat276', 'data_mat277', 'data_mat278', 'data_mat279', 'data_mat280', 'data_mat281', 'data_mat282', 'data_mat283', 'data_mat284', 'data_mat285', 'data_mat286', 'data_mat287', 'data_mat288', 'data_mat289', 'data_mat290', 'data_mat291', 'data_mat292', 'data_mat293', 'data_mat294', 'data_mat295', 'data_mat296', 'data_mat297', 'data_mat298', 'data_mat299', 'data_mat2100', 'data_mat2101', 'data_mat2102', 'data_mat2103', 'data_mat2104', 'data_mat2105', 'data_mat2106', 'data_mat2107', 'data_mat2108', 'data_mat2109', 'data_mat2110', 'data_mat2111', 'data_mat2112', 'data_mat2113', 'data_mat2114', 'data_mat2115', 'data_mat2116', 'data_mat2117', 'data_mat2118', 'data_mat2119', 'data_mat2120', 'data_mat2121', 'data_mat2122', 'data_mat2123', 'data_mat2124', 'data_mat2125', 'data_mat2126', 'data_mat2127', 'data_mat2128', 'data_mat2129', 'data_mat2130', 'data_mat2131', 'data_mat2132', 'data_mat2133', 'data_mat2134', 'data_mat2135', 'data_mat2136', 'data_mat2137', 'data_mat2138', 'data_mat2139', 'data_mat2140', 'data_mat2141', 'data_mat2142', 'data_mat2143', 'data_mat2144', 'data_mat2145', 'data_mat2146', 'data_mat2147', 'data_mat2148', 'data_mat2149', 'data_mat2150', 'data_mat2151', 'data_mat2152', 'data_mat2153', 'data_mat2154', 'data_mat2155', 'data_mat2156', 'data_mat2157', 'data_mat2158', 'data_mat2159', 'data_mat2160', 'data_mat2161', 'data_mat2162', 'data_mat2163', 'data_mat2164', 'data_mat2165', 'data_mat2166', 'data_mat2167', 'data_mat2168', 'data_mat2169', 'data_mat2170', 'data_mat2171', 'data_mat2172', 'data_mat2173', 'data_mat2174', 'data_mat2175', 'data_mat2176', 'data_mat2177', 'data_mat2178', 'data_mat2179', 'data_mat2180'};
trainedClassifier.ClassificationLinear = classificationLinear;
trainedClassifier.About = 'This struct is a trained model exported from Classification Learner R2023b.';
trainedClassifier.HowToPredict = sprintf('To make predictions on a new table, T, use: \n  [yfit,scores] = c.predictFcn(T) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nThe table, T, must contain the variables returned by: \n  c.RequiredVariables \nVariable formats (e.g. matrix/vector, datatype) must match the original training data. \nAdditional variables are ignored. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appclassification_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
inputTable = T_new;
predictorNames = {'data_mat21', 'data_mat22', 'data_mat23', 'data_mat24', 'data_mat25', 'data_mat26', 'data_mat27', 'data_mat28', 'data_mat29', 'data_mat210', 'data_mat211', 'data_mat212', 'data_mat213', 'data_mat214', 'data_mat215', 'data_mat216', 'data_mat217', 'data_mat218', 'data_mat219', 'data_mat220', 'data_mat221', 'data_mat222', 'data_mat223', 'data_mat224', 'data_mat225', 'data_mat226', 'data_mat227', 'data_mat228', 'data_mat229', 'data_mat230', 'data_mat231', 'data_mat232', 'data_mat233', 'data_mat234', 'data_mat235', 'data_mat236', 'data_mat237', 'data_mat238', 'data_mat239', 'data_mat240', 'data_mat241', 'data_mat242', 'data_mat243', 'data_mat244', 'data_mat245', 'data_mat246', 'data_mat247', 'data_mat248', 'data_mat249', 'data_mat250', 'data_mat251', 'data_mat252', 'data_mat253', 'data_mat254', 'data_mat255', 'data_mat256', 'data_mat257', 'data_mat258', 'data_mat259', 'data_mat260', 'data_mat261', 'data_mat262', 'data_mat263', 'data_mat264', 'data_mat265', 'data_mat266', 'data_mat267', 'data_mat268', 'data_mat269', 'data_mat270', 'data_mat271', 'data_mat272', 'data_mat273', 'data_mat274', 'data_mat275', 'data_mat276', 'data_mat277', 'data_mat278', 'data_mat279', 'data_mat280', 'data_mat281', 'data_mat282', 'data_mat283', 'data_mat284', 'data_mat285', 'data_mat286', 'data_mat287', 'data_mat288', 'data_mat289', 'data_mat290', 'data_mat291', 'data_mat292', 'data_mat293', 'data_mat294', 'data_mat295', 'data_mat296', 'data_mat297', 'data_mat298', 'data_mat299', 'data_mat2100', 'data_mat2101', 'data_mat2102', 'data_mat2103', 'data_mat2104', 'data_mat2105', 'data_mat2106', 'data_mat2107', 'data_mat2108', 'data_mat2109', 'data_mat2110', 'data_mat2111', 'data_mat2112', 'data_mat2113', 'data_mat2114', 'data_mat2115', 'data_mat2116', 'data_mat2117', 'data_mat2118', 'data_mat2119', 'data_mat2120', 'data_mat2121', 'data_mat2122', 'data_mat2123', 'data_mat2124', 'data_mat2125', 'data_mat2126', 'data_mat2127', 'data_mat2128', 'data_mat2129', 'data_mat2130', 'data_mat2131', 'data_mat2132', 'data_mat2133', 'data_mat2134', 'data_mat2135', 'data_mat2136', 'data_mat2137', 'data_mat2138', 'data_mat2139', 'data_mat2140', 'data_mat2141', 'data_mat2142', 'data_mat2143', 'data_mat2144', 'data_mat2145', 'data_mat2146', 'data_mat2147', 'data_mat2148', 'data_mat2149', 'data_mat2150', 'data_mat2151', 'data_mat2152', 'data_mat2153', 'data_mat2154', 'data_mat2155', 'data_mat2156', 'data_mat2157', 'data_mat2158', 'data_mat2159', 'data_mat2160', 'data_mat2161', 'data_mat2162', 'data_mat2163', 'data_mat2164', 'data_mat2165', 'data_mat2166', 'data_mat2167', 'data_mat2168', 'data_mat2169', 'data_mat2170', 'data_mat2171', 'data_mat2172', 'data_mat2173', 'data_mat2174', 'data_mat2175', 'data_mat2176', 'data_mat2177', 'data_mat2178', 'data_mat2179', 'data_mat2180'};
predictors = inputTable(:, predictorNames);
response = inputTable.Instrument;
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false];
classNames = [1; 2];

% Perform cross-validation
KFolds = 5;
cvp = cvpartition(response, 'KFold', KFolds);
% Initialize the predictions to the proper sizes
validationPredictions = response;
numObservations = size(predictors, 1);
numClasses = 2;
validationScores = NaN(numObservations, numClasses);
for fold = 1:KFolds
	trainingPredictors = predictors(cvp.training(fold), :);
	trainingResponse = response(cvp.training(fold), :);
	foldIsCategoricalPredictor = isCategoricalPredictor;

	% Train a classifier
	% This code specifies all the classifier options and trains the classifier.
	classificationLinear = fitclinear(...
		trainingPredictors, ...
		trainingResponse, ...
		'Learner', 'SVM', ...
		'Lambda', 'auto', ...
		'BetaTolerance', 0.0001, ...
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

%% Beta weights


%% Plot model evaluations 

% Compute classification using all data 
figure
C = confusionmat(response, validationPredictions);
confusionchart(C)

% Calculate accuracy
chart = confusionchart(response,validationPredictions); % Generate confusion chart
confusionMatrix = chart.NormalizedValues; % Get the normalized confusion matrix
accuracy(1) = sum(diag(confusionMatrix)) / sum(confusionMatrix(:)); % Calculate accuracy
title(sprintf('Accuracy = %0.2f%%', accuracy(1)*100))

% Plot all rates/VS 
figure
tiledlayout(2, 1)
mean_bass = mean(data_mat(1:20,:),1);
mean_oboe = mean(data_mat(21:40,:),1);
for ii = 1:2
	nexttile
	if ii == 1
		scatter(CFs/1000, mean_bass, 'blue')
		title('Bassoon')
	else
		scatter(CFs/1000, mean_oboe, 'red')
		title('Oboe')
	end
	set(gca, 'xscale', 'log')
	xticks([0.1 0.2 0.5 1 2 5 10])
	xlim([0.1 10])
	ylabel('Norm Rate')
	xlabel('CF (kHz)')
	grid on
end

%%

figure

nexttile
beta_weights = classificationLinear.Beta;
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
%xlim([0.2 10])
xlabel('CF')
ylabel('Beta Weight')
set(gca, 'xscale', 'log')

nexttile
scatter(CFs/1000, abs(beta_weights))
hold on 
yline(0)
xticks([0.1 0.2 0.5 1 2 5 10])
%xlim([0.2 10])
xlabel('CF')
ylabel('abs(Beta Weight)')
set(gca, 'xscale', 'log')

figure
tiledlayout(2, 1)
mean_bass = mean(data_mat(1:20,:),1);
mean_oboe = mean(data_mat(21:40,:),1);
for ii = 1:2
	nexttile
	if ii == 1
		scatter(CFs/1000, mean_bass, [], beta_weights, 'filled', ...
			'MarkerEdgeColor','k')
		title('Bassoon')
	else
		scatter(CFs/1000, mean_oboe, [], beta_weights, 'filled', ...
			'MarkerEdgeColor','k')
		title('Oboe')
	end
	set(gca, 'xscale', 'log')
	xticks([0.1 0.2 0.5 1 2 5 10])
	%xlim([0.2 10])
	ylabel('Norm Rate')
	xlabel('CF (kHz)')
	grid on
	colorbar
end

figure
tiledlayout(1, 2)
mean_bass = mean(data_mat(1:20,:),1);
mean_oboe = mean(data_mat(21:40,:),1);
index = abs(beta_weights)>0.4;
for ii = 1:2
	nexttile
	if ii == 1
		scatter(CFs(index)/1000, mean_bass(index), [], beta_weights(index), 'filled', ...
			'MarkerEdgeColor','k')
		title('Bassoon')
	else
		scatter(CFs(index)/1000, mean_oboe(index), [], beta_weights(index), 'filled', ...
			'MarkerEdgeColor','k')
		title('Oboe')
	end
	set(gca, 'xscale', 'log')
	xticks([0.1 0.2 0.5 1 2 5 10])
	%xlim([0.2 10])
	ylim([0 80])
	ylabel('Norm Rate')
	xlabel('CF (kHz)')
	grid on
	colorbar
end

%%

figure
tiledlayout(1, 2)
nexttile
scatter(mean_bass, mean_oboe, [], beta_weights, 'filled','MarkerEdgeColor','k')
hold on
plot([0 100], [0 100], 'k')
xlabel('Bassoon Norm Rate')
ylabel('Oboe Norm Rate')
xlim([0 100])
ylim([0 100])
colorbar

nexttile
scatter(mean_bass(index), mean_oboe(index), 'filled','MarkerEdgeColor','k')
hold on
plot([0 100], [0 100], 'k')
xlabel('Bassoon Norm Rate')
ylabel('Oboe Norm Rate')
xlim([0 100])
ylim([0 100])
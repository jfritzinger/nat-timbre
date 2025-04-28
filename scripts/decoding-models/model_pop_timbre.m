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

% Sort data 
for target = 1:16

	data_mat = NaN(2*20, num_data);
	for ii = 1:num_data

		if strcmp(type, 'rate')
			X1 = nat_data(sesh(ii)).bass_raterep(:,ind_b(target));
			X2 = nat_data(sesh(ii)).oboe_raterep(:,ind_o(target));
			X = [X1; X2];
		else
			X1 = nat_data(sesh(ii)).bass_VSrep(ind_b(target),:)';
			X2 = nat_data(sesh(ii)).oboe_VSrep(ind_o(target),:)';
			X = [X1; X2];
			X(X>=0.99) = NaN;
			X(isnan(X)) = 0;
		end
		data_mat(:,ii) = X;
	end

	% Create table for model
	T = array2table(data_mat);
	T.Response = response;
end

%% Run model with kfold validation 

% Take out test data
ncond = 2;
ind_test = false(ncond*20, 1); % Preallocate a logical array for 800 elements (40 * 20)
for istim = 1:ncond
	index = randperm(20, 4); % Randomly select 4 indices from 1 to 20
	ind_test((istim-1)*20 + index) = true; % Set the selected indices to true
end

% Make training and test rows
test_mat = data_mat(ind_test,:);
train_mat = data_mat(~ind_test,:);
test_response = [ones(1, 4) repmat(2, 1, 4)]';
train_response = [ones(1, 16) repmat(2, 1, 16)]';

T_train = array2table(train_mat);
T_train.Response = train_response;
T_test = array2table(test_mat);
T_test.Response = test_response;

%inputTable = T;
inputTable	= T_train;
predictors = inputTable(:, 1:180);
response = inputTable.Response;
classNames = [1; 2];

% Train a classifier
% This code specifies all the classifier options and trains the classifier.
classificationSVM = fitcsvm(...
	predictors, ...
	response, ...
	'KernelFunction', 'linear', ...
	'PolynomialOrder', [], ...
	'KernelScale', 'auto', ...
	'BoxConstraint', 1, ...
	'Standardize', true, ...
	'ClassNames', classNames);

% Create the result struct with predict function
predictorExtractionFcn = @(t) t(:,  1:180);
svmPredictFcn = @(x) predict(classificationSVM, x);
trainedClassifier.predictFcn = @(x) svmPredictFcn(predictorExtractionFcn(x));

% Add additional fields to the result struct
trainedClassifier.RequiredVariables = {'data_mat1', 'data_mat2', 'data_mat3', 'data_mat4', 'data_mat5', 'data_mat6', 'data_mat7', 'data_mat8', 'data_mat9', 'data_mat10', 'data_mat11', 'data_mat12', 'data_mat13', 'data_mat14', 'data_mat15', 'data_mat16', 'data_mat17', 'data_mat18', 'data_mat19', 'data_mat20', 'data_mat21', 'data_mat22', 'data_mat23', 'data_mat24', 'data_mat25', 'data_mat26', 'data_mat27', 'data_mat28', 'data_mat29', 'data_mat30', 'data_mat31', 'data_mat32', 'data_mat33', 'data_mat34', 'data_mat35', 'data_mat36', 'data_mat37', 'data_mat38', 'data_mat39', 'data_mat40', 'data_mat41', 'data_mat42', 'data_mat43', 'data_mat44', 'data_mat45', 'data_mat46', 'data_mat47', 'data_mat48', 'data_mat49', 'data_mat50', 'data_mat51', 'data_mat52', 'data_mat53', 'data_mat54', 'data_mat55', 'data_mat56', 'data_mat57', 'data_mat58', 'data_mat59', 'data_mat60', 'data_mat61', 'data_mat62', 'data_mat63', 'data_mat64', 'data_mat65', 'data_mat66', 'data_mat67', 'data_mat68', 'data_mat69', 'data_mat70', 'data_mat71', 'data_mat72', 'data_mat73', 'data_mat74', 'data_mat75', 'data_mat76', 'data_mat77', 'data_mat78', 'data_mat79', 'data_mat80', 'data_mat81', 'data_mat82', 'data_mat83', 'data_mat84', 'data_mat85', 'data_mat86', 'data_mat87', 'data_mat88', 'data_mat89', 'data_mat90', 'data_mat91', 'data_mat92', 'data_mat93', 'data_mat94', 'data_mat95', 'data_mat96', 'data_mat97', 'data_mat98', 'data_mat99', 'data_mat100', 'data_mat101', 'data_mat102', 'data_mat103', 'data_mat104', 'data_mat105', 'data_mat106', 'data_mat107', 'data_mat108', 'data_mat109', 'data_mat110', 'data_mat111', 'data_mat112', 'data_mat113', 'data_mat114', 'data_mat115', 'data_mat116', 'data_mat117', 'data_mat118', 'data_mat119', 'data_mat120', 'data_mat121', 'data_mat122', 'data_mat123', 'data_mat124', 'data_mat125', 'data_mat126', 'data_mat127', 'data_mat128', 'data_mat129', 'data_mat130', 'data_mat131', 'data_mat132', 'data_mat133', 'data_mat134', 'data_mat135', 'data_mat136', 'data_mat137', 'data_mat138', 'data_mat139', 'data_mat140', 'data_mat141', 'data_mat142', 'data_mat143', 'data_mat144', 'data_mat145', 'data_mat146', 'data_mat147', 'data_mat148', 'data_mat149', 'data_mat150', 'data_mat151', 'data_mat152', 'data_mat153', 'data_mat154', 'data_mat155', 'data_mat156', 'data_mat157', 'data_mat158', 'data_mat159', 'data_mat160', 'data_mat161', 'data_mat162', 'data_mat163', 'data_mat164', 'data_mat165', 'data_mat166', 'data_mat167', 'data_mat168', 'data_mat169', 'data_mat170', 'data_mat171', 'data_mat172', 'data_mat173', 'data_mat174', 'data_mat175', 'data_mat176', 'data_mat177', 'data_mat178', 'data_mat179', 'data_mat180'};
trainedClassifier.ClassificationSVM = classificationSVM;
trainedClassifier.About = ['This struct is a trained model exported from ' ...
	'Classification Learner R2023b.'];
trainedClassifier.HowToPredict = sprintf(['To make predictions on a new table,' ...
	' T, use: \n  [yfit,scores] = c.predictFcn(T) \nreplacing ''c'' with ' ...
	'the name of the variable that is this struct, e.g. ''trainedModel''.' ...
	' \n \nThe table, T, must contain the variables returned by: \n  ' ...
	'c.RequiredVariables \nVariable formats (e.g. matrix/vector, datatype) ' ...
	'must match the original training data. \nAdditional variables are ignored.' ...
	' \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot,' ...
	' ''stats'', ''stats.map''), ''appclassification_exportmodeltoworkspace'')"' ...
	'>How to predict using an exported model</a>.']);

% Perform cross-validation
partitionedModel = crossval(trainedClassifier.ClassificationSVM, 'KFold', 5);

% Compute validation predictions
[validationPredictions, validationScores] = kfoldPredict(partitionedModel);

% Compute validation accuracy
validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');

% % This code processes the data into the right shape for testing the
% % model.
% inputTable = T_test;
% predictors = inputTable(:, 1:180);
% response = inputTable.Response;
% classNames = [1; 2];
% 
% % Run model with testing data 
% [yfit,scores] = trainedClassifier.predictFcn(inputTable);


%% Plot model evaluations 

% Compute classification using all data 
figure
C = confusionmat(train_response, validationPredictions);
confusionchart(C)

% Calculate accuracy
chart = confusionchart(train_response,validationPredictions); % Generate confusion chart
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
beta_weights = classificationSVM.Beta;
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
	xlim([0.2 10])
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
	xlim([0.2 10])
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

%% Run model
%
% nexttile
% hold on
% swarmchart(ones(10,1), alloutput_avg(:,1), 'filled')
% swarmchart(ones(10,1)*2, alloutput_avg(:,2), 'filled')
% xlim([0.5 2.5])
% ylim([0.5 2.5])
% xlabel('Actual Timbre')
% ylabel('Predicted Timbre')
%
% % Plot accuracy
% nexttile
% swarmchart(ones(10,1), accuracy_one(:,1), 'filled')
% hold on
% swarmchart(ones(10,1)*2, accuracy_one(:,2), 'filled')
% xlim([0.5 2.5])
% %set(gca, 'XScale', 'log')

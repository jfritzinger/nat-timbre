%% model_neuron_rate_F0
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
[F0s1, order] = sort(F0s1);
F0s1 = log10(F0s1);
response = reshape(repmat(F0s1, 1, 20)', 1, []);
response = response';

%% Get all rates for each repetition for bassoon (one example neuron)

%putative = 'R29_TT2_P3_N03'; % Predictions 44, yesterday 62? 
% putative = 'R29_TT3_P5_N05'; % Not predicting well, not sure why 
putative = 'R29_TT3_P2_N06'; % Good 
index = find(strcmp({nat_data.putative}, putative));

%% Run model 
rng("shuffle"); % Seed based on current time

figure
tiledlayout(5, 2)
for imodelrep = 1:10

	% Get data
	h_all = [];
	for itarget = 1:40
		spikes_bass = nat_data(index).bass_spikerate{itarget}/1000; % ms
		spikereps_bass = nat_data(index).bass_spikerep{itarget};

		% Arrange data for SVM
		min_dis = 0.25;
		edges = 0:min_dis:300;
		t = 0+min_dis/2:min_dis:300-min_dis/2;
		for irep = 1:20
			%h_bass(irep, :) = histcounts(spikes_bass(spikereps_bass==irep), edges);
			rep_psth = histcounts(spikes_bass(spikereps_bass==irep), edges);
			shuffled_data = rep_psth(randperm(numel(rep_psth)));
			h_bass(irep, :) = shuffled_data;
		end
		h_all = [h_all; h_bass];
	end

	% Put data into table
	T = array2table(h_all);
	T.response = response;
	predictors = h_all;

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
	numClasses = 40;
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
		validationPredictions (cvp.test(fold), :) = foldPredictions;
		validationScores(cvp.test(fold), :) = foldScores;
	end

	% Compute validation accuracy
	correctPredictions = (validationPredictions == response);
	isMissing = isnan(response);
	correctPredictions = correctPredictions(~isMissing);
	accur(imodelrep) = sum(correctPredictions)/length(correctPredictions);
	C = confusionmat(validationPredictions, response);

	% Get accuracy for 1-20, 21-40
	accur1 = sum(correctPredictions(1:400))/length(correctPredictions(1:400));
	accur2 = sum(correctPredictions(401:800))/length(correctPredictions(401:800));

	nexttile
	confusionchart(C)
	title(num2str(accur(imodelrep)*100))
end

%% Run model (normal)

figure
for imodelrep = 1

	% Get data
	h_all = [];
	for itarget = 1:40
		spikes_bass = nat_data(index).bass_spikerate{itarget}/1000; % ms
		spikereps_bass = nat_data(index).bass_spikerep{itarget};

		% Arrange data for SVM
		min_dis = 0.25;
		edges = 0:min_dis:300;
		t = 0+min_dis/2:min_dis:300-min_dis/2;
		for irep = 1:20
			h_bass(irep, :) = histcounts(spikes_bass(spikereps_bass==irep), edges);
		end
		h_all = [h_all; h_bass];
	end

	% Put data into table
	T = array2table(h_all);
	T.response = response;
	predictors = h_all;

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
	numClasses = 40;
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

	% Compute validation accuracy
	correctPredictions = (validationPredictions == response);
	isMissing = isnan(response);
	correctPredictions = correctPredictions(~isMissing);
	accur(imodelrep) = sum(correctPredictions)/length(correctPredictions);
	C = confusionmat(validationPredictions, response);

	% Get accuracy for 1-20, 21-40
	accur1 = sum(correctPredictions(1:400))/length(correctPredictions(1:400));
	accur2 = sum(correctPredictions(401:800))/length(correctPredictions(401:800));

	nexttile
	confusionchart(C)
	title(num2str(accur(imodelrep)*100))
end

%% Plot predicted vs actual 

figure('Position',[38,937,880,332])
tiledlayout(1, 2)
nexttile
scatter(response, validationPredictions, 'filled', 'MarkerFaceAlpha',0.2)
set(gca, 'xscale', 'log', 'yscale', 'log')
ylim([57 588])
xlim([57 588])
xticks([58 110 220 440 580])
yticks([58 110 220 440 580])
grid on
title(['Accuracy = ' num2str(accur(imodelrep)*100) '%'])


nexttile
for ii = 1:40
	start_idx = (ii-1)*20 + 1;
    end_idx = ii*20;
	mean_resp(ii) = mean(response(start_idx:end_idx));
	mean_pred(ii) = mean(validationPredictions(start_idx:end_idx));
end
plot(mean_resp, mean_resp)
hold on
scatter(mean_resp, mean_pred, 'filled')
set(gca, 'xscale', 'log', 'yscale', 'log')
ylim([57 588])
xlim([57 588])
r = corrcoef(mean_resp, mean_pred);
r2 = r(1, 2)^2;
title(['R^2 = ' num2str(r2)])
xticks([58 110 220 440 580])
yticks([58 110 220 440 580])
box off
grid on

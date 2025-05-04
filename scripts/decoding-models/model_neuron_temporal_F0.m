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

F0s = log10(F0s1);
response = reshape(repmat(F0s, 1, 20)', 1, []);
response = response';

%% Get all rates for each repetition for bassoon (one example neuron)

sesh = find(~cellfun(@isempty, {nat_data.bass_rate}));
num_data = numel(sesh);

for ind = 1:num_data
	index = sesh(ind);

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
	for imodelrep = 1 %:5

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
		validationAccuracy(ind) = sum(correctPredictions)/length(correctPredictions);
		C = confusionmat(validationPredictions, response);

		% Get accuracy for 1-20, 21-40
		validationAccuracy1(ind) = sum(correctPredictions(1:400))/length(correctPredictions(1:400));
		validationAccuracy2(ind) = sum(correctPredictions(401:800))/length(correctPredictions(401:800));
		
	end

	figure
	confusionchart(C)
	title(num2str(validationAccuracy*100))
	
	% Save data for each
	neuron_time_F0(index).putative = nat_data(index).putative;
	neuron_time_F0(index).CF = nat_data(index).CF;
	neuron_time_F0(index).MTF = nat_data(index).MTF;
	neuron_time_F0(index).response = response;
	neuron_time_F0(index).predictors = predictors;
	neuron_time_F0(index).T = T;
	neuron_time_F0(index).validationPredictions = validationPredictions;
	neuron_time_F0(index).accuracy = validationAccuracy;
	neuron_time_F0(index).accuracy_low = validationAccuracy1;
	neuron_time_F0(index).accuracy_high = validationAccuracy2;
	neuron_time_F0(index).C = C;
	
	fprintf('%d/%d, %0.2f%% done!\n', ind, num_data, ind/num_data*100)
end

%% Save struct of data 

save('Neuron_Time_F0_Bassoon.mat', "neuron_time_F0")


%% Plot accuracy of each neuron
figure('Position',[560,618,798,230])
tiledlayout(1, 3)

nexttile
edges = linspace(0, 100, 30);
histogram(validationAccuracy*100,edges)
mean_F0 = mean(validationAccuracy);
hold on
xline(2.5, 'k')
xline(mean_F0*100, 'r', 'LineWidth',2)
ylabel('# Neurons')
xlabel('Prediction Accuracy (%)')
title('Prediction of F0')
mean_all = mean(validationAccuracy, 'all');
fprintf('Mean for all = %0.4f\n', mean_all)
xlim([0 100])
legend('', 'Chance', ['Mean=' num2str(round(mean_all*100)) '%'])

nexttile
histogram(validationAccuracy1*100,edges)
mean_F0 = mean(validationAccuracy1);
hold on
xline(2.5, 'k')
xline(mean_F0*100, 'r', 'LineWidth',2)
ylabel('# Neurons')
xlabel('Prediction Accuracy (%)')
title('Prediction of F0, 58-174Hz')
mean_all = mean(validationAccuracy1, 'all');
fprintf('Mean for all = %0.4f\n', mean_all)
xlim([0 100])
legend('', 'Chance', ['Mean=' num2str(round(mean_all*100)) '%'])

nexttile
histogram(validationAccuracy2*100,edges)
mean_F0 = mean(validationAccuracy2);
hold on
xline(2.5, 'k')
xline(mean_F0*100, 'r', 'LineWidth',2)
ylabel('# Neurons')
xlabel('Prediction Accuracy (%)')
title('Prediction of F0, 185-587Hz')
mean_all = mean(validationAccuracy2, 'all');
fprintf('Mean for all = %0.4f\n', mean_all)
xlim([0 100])
legend('', 'Chance', ['Mean=' num2str(round(mean_all*100)) '%'])


%% Find indices of the best neurons 

[temp,originalpos] = sort(validationAccuracy, 'descend' );
n = temp(1:3);
best_ind=originalpos(1:3);

putatives = nat_data(best_ind(1)).putative;

%%

for ind = best_ind
	index = sesh(ind);

	% Get data
	h_all = [];
	for itarget = 1:40
		spikes_bass = nat_data(index).bass_spikerate{itarget}/1000; % ms
		spikereps_bass = nat_data(index).bass_spikerep{itarget};

		% Arrange data for SVM
		min_dis = 0.5;
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
	for imodelrep = 1 %:5

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
		accur(ind) = sum(correctPredictions)/length(correctPredictions);
		C = confusionmat(validationPredictions, response);

		% Get accuracy for 1-20, 21-40
		accur1(ind) = sum(correctPredictions(1:400))/length(correctPredictions(1:400));
		accur2(ind) = sum(correctPredictions(401:800))/length(correctPredictions(401:800));
		
	end

	figure
	confusionchart(C)
	title(num2str(accur(ind)*100))
end

%% Plot rasters of three best




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
	T.response = response';

	predictors = h_all;
	response = response';
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
		validationAccuracy = sum(correctPredictions)/length(correctPredictions);

		C = confusionmat(validationPredictions, response);
		
	end

	figure
	confusionchart(C)
	title(num2str(validationAccuracy*100))


	% Plots
	% figure
	% tiledlayout(2, 1)
	% nexttile
	% histogram(spikes_bass, 301)
	% nexttile
	% histogram(spikes_oboe, 301)
	% figure
	% tiledlayout(2, 1)
	% nexttile
	% hold on
	% for irep = 1:20
	% 	plot(t,h_bass(irep,:)+irep)
	% end
	% scatter(spikes_bass, spikereps_bass, 'filled')
	% nexttile
	% hold on
	% for irep = 1:20
	% 	plot(t,h_oboe(irep,:)+irep)
	% end
	% scatter(spikes_oboe, spikereps_oboe, 'filled')
	% 
	% fprintf('%d/%d, %0.2f%% done!\n', ind, num_data, ind/num_data*100)
end

% Plot accuracy of each neuron
figure
tiledlayout(4, 4)
acc2 = squeeze(max(accuracy, 3));
for ii = 1:16
	nexttile
	histogram(acc2(:,ii)*100,21)
	mean_F0 = mean(acc2(:,ii));
	hold on
	xline(50, 'k')
	xline(mean_F0*100, 'r', 'LineWidth',2)
	ylabel('# Neurons')
	xlabel('Prediction Accuracy (%)')
	title(['Prediction of instrument, F0=' num2str(round(bass_pitch(ind_b(ii))))])
end

mean_all = mean(acc2, 'all');
fprintf('Mean for all = %0.4f\n', mean_all)


%% model_neuron_rate_F0
clear

%% Load in data

if ismac
	filepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/model_comparisons';
else
	filepath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\model_comparisons';
end
load(fullfile(filepath, 'Data_NT.mat'), 'nat_data')


%% Get correct output of model 
target = 'Oboe';

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
[F0s1, order] = sort(F0s1);

F0s = log10(F0s1);
response = reshape(repmat(F0s, 1, 20)', 1, []);
response = response';

%% Get all rates for each repetition for bassoon (one example neuron)

sesh = find(~cellfun(@isempty, {nat_data.oboe_rate}));
num_data = numel(sesh);

for ind = 1:num_data
	index = sesh(ind);

	% Get data
	h_all = [];
	for itarget = 1:length(F0s)
		if strcmp(target, 'Oboe')
			spikes = nat_data(index).oboe_spikerate{itarget}/1000; % ms
			spikereps = nat_data(index).oboe_spikerep{itarget};
		else
			spikes = nat_data(index).bass_spikerate{itarget}/1000; % ms
			spikereps = nat_data(index).bass_spikerep{itarget};
		end

		% Arrange data for SVM
		min_dis = 0.25;
		edges = 0:min_dis:300;
		t = 0+min_dis/2:min_dis:300-min_dis/2;
		for irep = 1:20
			h_bass(irep, :) = histcounts(spikes(spikereps==irep), edges);
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

		% Compute validation accuracy
		correctPredictions = (validationPredictions == response);
		isMissing = isnan(response);
		correctPredictions = correctPredictions(~isMissing);
		validationAccuracy(ind) = sum(correctPredictions)/length(correctPredictions);
		C = confusionmat(validationPredictions, response);

		% Get accuracy for 1-20, 21-40
		validationAccuracy1(ind) = sum(correctPredictions(1:400))/length(correctPredictions(1:400));
		validationAccuracy2(ind) = sum(correctPredictions(401:end))/length(correctPredictions(401:end));
		
	end

% 	figure
% 	confusionchart(C)
% 	title(num2str(validationAccuracy*100))
	
	% Save data for each
	neuron_time_F0(ind).putative = nat_data(index).putative;
	neuron_time_F0(ind).CF = nat_data(index).CF;
	neuron_time_F0(ind).MTF = nat_data(index).MTF;
	neuron_time_F0(ind).response = response;
	neuron_time_F0(ind).predictors = predictors;
	neuron_time_F0(ind).T = T;
	neuron_time_F0(ind).validationPredictions = validationPredictions;
	neuron_time_F0(ind).accuracy = validationAccuracy(ind);
	neuron_time_F0(ind).accuracy_low = validationAccuracy1(ind);
	neuron_time_F0(ind).accuracy_high = validationAccuracy2(ind);
	neuron_time_F0(ind).C = C;
	
	fprintf('%d/%d, %0.2f%% done!\n', ind, num_data, ind/num_data*100)
end

%% Save struct of data 

[base, datapath, savepath, ppi] = getPathsNT();
save(fullfile(base, 'model_comparisons', ['Neuron_Time_F0_' target '.mat']), ...
	"neuron_time_F0", '-v7.3')

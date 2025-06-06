%% model_neuron_rate_F0
clear

%% Load in data

[base, datapath, savepath, ppi] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')


%% Get correct output of model 
target = 'Bassoon';

% Get bassoon stimulus
tuning = readtable(fullfile(base, 'Tuning.xlsx')); % Load in tuning
listing = dir(fullfile(base, 'waveforms', ['*' target '*.wav']));
files = {listing.name};
note_names = extractBetween(files, 'ff.', '.');
[~, index] = ismember(note_names, tuning.Note);
F0s1 = tuning.Frequency(index);
[F0s1, order] = sort(F0s1);

F0s = log10(F0s1);
response = reshape(repmat(F0s, 1, 20)', 1, []);
response = response';

%% Get all rates for each repetition for bassoon (one example neuron)

if strcmp(target, 'Oboe')
	sesh = find(~cellfun(@isempty, {nat_data.oboe_rate}));
else
	sesh = find(~cellfun(@isempty, {nat_data.bass_rate}));
end
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

	for imodelrep = 1 %:5

		[validationPredictions] = trainClassifierNeuronTimeF0(T);

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

save(fullfile(base, 'model_comparisons', ['Neuron_Time_F0_' target '.mat']), ...
	"neuron_time_F0", '-v7.3')

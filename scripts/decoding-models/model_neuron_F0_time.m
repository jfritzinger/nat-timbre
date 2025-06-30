%% model_neuron_rate_F0
clear

%% Load in data

[base, datapath, savepath, ppi] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')


%% Get correct output of model 
%target = 'Bassoon';
%target = 'Oboe';
target = 'Invariant';

% Get bassoon stimulus
F0s = getF0s(target);
F0s = log10(F0s);

%% Get all rates for each repetition for bassoon (one example neuron)

[sesh, num_data] = getF0Sessions(nat_data, target);
for ind = 1:num_data
	index = sesh(ind);
	T = getF0NeuronTable(nat_data, target, index, F0s, 'Timing');

	for imodelrep = 1 %:5

		[validationPredictions] = trainClassifierNeuronTimeF0(T, F0s);

		% Compute validation accuracy
		correctPredictions = (validationPredictions == T.Response);
		validationAccuracy(ind) = sum(correctPredictions)/length(correctPredictions);
		C = confusionmat(validationPredictions, T.Response);
	end

% 	figure
% 	confusionchart(C)
% 	title(num2str(validationAccuracy*100))
	
	% Save data for each
	neuron_time_F0(ind).putative = nat_data(index).putative;
	neuron_time_F0(ind).CF = nat_data(index).CF;
	neuron_time_F0(ind).MTF = nat_data(index).MTF;
	neuron_time_F0(ind).response = T.Response;
	neuron_time_F0(ind).predictors = predictors;
	neuron_time_F0(ind).T = T;
	neuron_time_F0(ind).validationPredictions = validationPredictions;
	neuron_time_F0(ind).accuracy = validationAccuracy(ind);
	neuron_time_F0(ind).C = C;
	
	fprintf('%d/%d, %0.2f%% done!\n', ind, num_data, ind/num_data*100)
end

%% Save struct of data 

save(fullfile(base, 'model_comparisons', ['Neuron_Time_F0_' target '.mat']), ...
	"neuron_time_F0", '-v7.3')

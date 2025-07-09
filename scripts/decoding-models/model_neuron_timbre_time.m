%% model_neuron_rate_F0
clear


%% Load in data

[base, ~, ~, ~] = getPathsNT();
%load(fullfile(base,'model_comparisons', 'Data_NT_3.mat'), 'nat_data')
load(fullfile(base, 'model_comparisons',  'Model_NT2.mat'), 'nat_model')
nat_data = nat_model;

%% Shape data into model input

[sesh, num_data] = getTimbreSessions(nat_data);
for ind = 1:num_data
	index = sesh(ind);

	T = getTimbreNeuronTable(nat_data, index, 'Model');

	% Call SVM
	[trainedClassifier, validationAccuracy, validationPredictions] = ...
		trainClassifierTimeTimbre(T);
	C = confusionmat(validationPredictions, T.Instrument);
	accuracy_all(ind) = sum(diag(C)) / sum(C, "all");

	% Set up struct to save data
	neuron_time_timbre(ind).putative = nat_data(index).putative;
	neuron_time_timbre(ind).CF = nat_data(index).CF;
	neuron_time_timbre(ind).MTF = nat_data(index).MTF;
	neuron_time_timbre(ind).T = T;
	neuron_time_timbre(ind).Instrument = T.Instrument;
	neuron_time_timbre(ind).prediction = validationPredictions;
	neuron_time_timbre(ind).accuracy = accuracy_all(ind);
	neuron_time_timbre(ind).C = C;

	fprintf('%d/%d, %0.2f%% done! Accur=%0.2f%%\n', ind, ...
		num_data, ind/num_data*100, accuracy_all(ind)*100)
end

% save(fullfile(base, 'model_comparisons', 'Neuron_Time_Timbre_All_Coarse300.mat'), ...
% 	"neuron_time_timbre", '-v7.3')
save(fullfile(base, 'model_comparisons', 'Model_N_Time_Timbre_All.mat'), ...
	"neuron_time_timbre", '-v7.3')

%% Plot accuracy of each neuron

figure
histogram(accuracy_all*100,51)
mean_F0 = mean(accuracy_all);
hold on
xline(50, 'k')
xline(mean_F0*100, 'r', 'LineWidth',2)
ylabel('# Neurons')
xlabel('Prediction Accuracy (%)')
xlim([0 100])
mean_all = mean(accuracy_all, 'all');
fprintf('Mean for all = %0.4f\n', mean_all)
title('Instrument identification task using SFIE BE/BS timing')
%save('Neuron_Time_Timbre_All.mat', "accuracy", "mean_F0")

%% Plot scatter plot 

[base, datapath, savepath, ppi] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Model_Neuron_Rate_Timbre_All.mat'), ...
	"neuron_rate_timbre")
accuracy_rate = [neuron_rate_timbre.accuracy];

figure
scatter(accuracy_rate, accuracy_all, 'filled', 'MarkerEdgeColor','k');
hold on
xlim([0.4 1])
ylim([0.4 1])
plot([0.4 1], [0.4 1], 'k')
ylabel('')
xlabel('Rate Accuracy')
ylabel('Timing Accuracy')
title('Model using shuffled              PSTH')

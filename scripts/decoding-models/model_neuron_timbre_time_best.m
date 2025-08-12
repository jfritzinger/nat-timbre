%% model_neuron_rate_F0
clear


%% Load in data

[base, ~, ~, ~] = getPathsNT();
load(fullfile(base,'model_comparisons_revised', 'Data_NT_3.mat'), 'nat_data')

%% Shape data into model input

[sesh, num_data] = getTimbreSessions(nat_data);
neuron_time_timbre = struct();

min_dis = [0.1 0.25 0.5 0.75 1 1.5 2 2.5 3 3.5 4 4.5 5 10 15 20 30 40 50 100 150];
accuracy_all = zeros(length(min_dis), num_data);
for idist = 1:length(min_dis)
	min_dis2 = min_dis(idist);
	for ind = 1:num_data
		index = sesh(ind);

		% Separate out training and testing data
		T = getTimbreNeuronTable(nat_data, index, 'Data', min_dis2);
		[T_train, T_test] = splitData(T);

		% Call SVM
		[trainedClassifier, validationAccuracy, validationPredictions] = ...
			trainClassifierTimeTimbre(T_train);

		% To make predictions with the returned 'trainedClassifier' on new data
		[yfit,scores] = trainedClassifier.predictFcn(T_test);
		C = confusionmat(T_test.Response, yfit);
		accuracy_all(idist, ind) = sum(diag(C))/sum(C, 'all');

		% Set up struct to save data
		neuron_time_timbre(idist, ind).putative = nat_data(index).putative;
		neuron_time_timbre(idist, ind).accuracy = accuracy_all(idist, ind);
		neuron_time_timbre(idist, ind).C = C;
		neuron_time_timbre(idist, ind).min_dis = min_dis2;

		fprintf('%d/%d, %0.2f%% done! Accur=%0.2f%%\n', ind, ...
			num_data, ind/num_data*100, accuracy_all(idist, ind)*100)
	end
end
save(fullfile(base, 'model_comparisons_revised', 'Neuron_Time_Timbre_Distances.mat'), ...
	"neuron_time_timbre", '-v7.3')
% save(fullfile(base, 'model_comparisons', 'Model_N_Time_Timbre_All.mat'), ...
% 	"neuron_time_timbre", '-v7.3')

%% Plot accuracy of each neuron

figure
tiledlayout(7, 3, 'TileSpacing','tight', 'Padding','compact')
for i = 1:length(min_dis)
	nexttile
	histogram(accuracy_all(i,:)*100,51)
	mean_F0 = mean(accuracy_all(i,:));
	hold on
	xline(50, 'k')
	xline(mean_F0*100, 'r', 'LineWidth',1)
	ylabel('# Neurons')
	xlabel('Prediction Accuracy (%)')
	xlim([0 100])
	mean_all = mean(accuracy_all(i,:), 'all');
	fprintf('Mean for all = %0.4f\n', mean_all)
	title(['Min Dist = ' num2str(min_dis(i)) ', Mean = ' num2str(round(mean_all*100))])
end
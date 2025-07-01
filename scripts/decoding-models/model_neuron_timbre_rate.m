%% model_neuron_rate_F0
clear

%% Load in data

[base, datapath, savepath, ppi] = getPathsNT();
%load(fullfile(filepath, 'Data_NT_3.mat'), 'nat_data')
load(fullfile(base, 'model_comparisons',  'Model_NT.mat'), 'nat_model')
nat_data = nat_model;

%% Shape data into model input

[sesh, num_data] = getTimbreSessions(nat_data);
ind_b = 25:40;
ind_o = [1 3:17];
target = 1;

% Get all rates for each repetition for bassoon (one example neuron)
for ind = 1:num_data

	index = sesh(ind);
	for target = 1:16
		data_1 = [nat_data(index).bass_raterep(:, ind_b(target)) ...
			nat_data(index).oboe_raterep(:, ind_o(target))];
		idx = (1:20) + 20*(target-1);
		data(idx, :) = data_1;
	end

	% Rearrange data 
	data_all = [data(:,1); data(:,2)];
	T = array2table(data_all);
	response = [ones(1, 320) repmat(2, 1, 320)]';
	T.Instrument = response;

	%% Calculate simple rate prediction

	%[trainedClassifier, accuracy_SVM(ind)] = trainClassifierNeuronTimbre(T);

	% Initialize variables to store results
	closest = zeros(size(data, 1), size(data, 2));
	actual = zeros(size(data, 1), size(data, 2));

	% Loop through each row to calculate the average of the other rows and find the closest column
	for i = 1:size(data, 1)

		% Extract the current row (1x40 data)
		single_row = data(i, :);

		% Calculate overall average rates across repetitions in response to each of
		% the 12 vowels, excluding the current repetition.
		other_rows = data([1:i-1, i+1:end], :); % Remove row i
		avg_other_rows = mean(other_rows, 1); % Average of other rows (1x40)

		% Calculate average rate to a given repetition of one vowel
		for ii = 1:2
			rate = single_row(ii);

			% The response to each repetition was identified as the vowel for which the
			% absolute difference between the single-repetition rate and overall
			% average rate was minimal.
			differences = abs(avg_other_rows - rate);
			[~, closest_column_index] = min(differences);

			% Store the closest column index for this row
			closest(i, ii) = closest_column_index;
			actual(i,ii) = ii;
		end
	end

	%% Analysis

	% Calculate accuracy
	actual2 = reshape(actual,[], 320*2);
	closest2 = reshape(closest, [], 320*2);
	C = confusionmat(actual2, closest2);
	chart = confusionchart(actual2,closest2); % Generate confusion chart
	confusionMatrix = chart.NormalizedValues; % Get the normalized confusion matrix
	accuracy(ind) = sum(diag(confusionMatrix)) / sum(confusionMatrix(:)); % Calculate accuracy

	% Plot average rates
	% nexttile
	% hold on
	% bar(avg_rate)
	% errorbar(1:2, avg_rate, rate_std/sqrt(20), 'LineStyle','none', 'Color','k')
	% ylabel('Avg. Rate')
	% xlabel('F0s')
	% title(sprintf('F0 = %0.0f, Acc = %0.1f%%', ...
	% 	bass_pitch(ind_b(target)), accuracy(ind, target)*100))
	%
	% Plot confusion matrix
	%nexttile
	%confusionchart(C)

	% Set up struct to save data
	neuron_rate_timbre(ind).putative = nat_data(index).putative;
	neuron_rate_timbre(ind).ind_b = ind_b(target);
	neuron_rate_timbre(ind).ind_o = ind_o(target);
	neuron_rate_timbre(ind).CF = nat_data(index).CF;
	neuron_rate_timbre(ind).MTF = nat_data(index).MTF;
	neuron_rate_timbre(ind).rate_rep = data;
	neuron_rate_timbre(ind).actual = actual2;
	neuron_rate_timbre(ind).closest = closest2;
	neuron_rate_timbre(ind).accuracy = accuracy(ind);
	neuron_rate_timbre(ind).C = C;

	fprintf('%d/%d, %0.2f%% done!\n', ind, num_data, ind/num_data*100)
end

%% Plot outputs 

edges = linspace(0, 1, 51);
figure
histogram(accuracy, edges)
ylabel('Number of neurons')
title('Instrument identification task using SFIE BE/BS rate')
xlabel('Accuracy')

% figure
% nexttile
% scatter(accuracy, accuracy_SVM, 'filled', 'MarkerFaceAlpha',0.5)
% hold on
% plot([0.4 1], [0.4 1])
% ylabel('SVM Accuracy')
% xlabel('Rate Discrimination Accuracy')
% title('Comparing SVM to manual discrimination')
% 
% % Load in timing models
% filepath_timing = fullfile(base, 'model_comparisons', 'Neuron_Time_Timbre_All.mat');
% load(filepath_timing, "neuron_time_timbre")
% accuracy_time = [neuron_time_timbre.accuracy];
% nexttile
% scatter(accuracy, accuracy_time, 'filled', 'MarkerFaceAlpha',0.5)
% hold on
% plot([0.4 1], [0.4 1])
% ylabel('Time')
% xlabel('Rate Discrimination Accuracy')
% 
% nexttile
% scatter(accuracy_SVM, accuracy_time, 'filled', 'MarkerFaceAlpha',0.5)
% hold on
% plot([0.4 1], [0.4 1])
% ylabel('Time')
% xlabel('SVM Rate Accuracy')

%% Save data

% save(fullfile(base, 'model_comparisons', 'Neuron_Rate_Timbre_All.mat'), ...
% 	"neuron_rate_timbre")
save(fullfile(base, 'model_comparisons', 'Model_N_Rate_Timbre_All.mat'), ...
	"neuron_rate_timbre")

%% model_neuron_rate_F0
clear 
%% Load in data

filepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/model_comparisons';
load(fullfile(filepath, 'Data_NT.mat'), 'nat_data')


%% Shape data into model input

sesh = find(~cellfun(@isempty, {nat_data.oboe_rate}));
num_data = numel(sesh);

% Get all rates for each repetition for bassoon (one example neuron)
neuron_rate_F0 = struct;
for ind = 1:num_data

	index = sesh(ind);
	data = nat_data(index).oboe_raterep;
	avg_rate = nat_data(index).oboe_rate;
	rate_std = nat_data(index).oboe_rate_std;

	%% Calculate simple rate prediction

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
		for ii = 1:length(avg_rate)
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

	% Plot average rates
	% figure('Position',[136,782,1085,481])
	% tiledlayout(1, 2);
	% nexttile
	% hold on
	% bar(avg_rate)
	% errorbar(1:40, avg_rate, rate_std/sqrt(20), 'LineStyle','none', 'Color','k')
	% ylabel('Avg. Rate')
	% xlabel('F0s')

	% Plot confusion matrix
	actual2 = reshape(actual,[], 20*length(avg_rate));
	closest2 = reshape(closest, [], 20*length(avg_rate));
	% nexttile
	C = confusionmat(actual2, closest2);
	%confusionchart(C)

	% Calculate accuracy
	chart = confusionchart(actual2,closest2); % Generate confusion chart
	confusionMatrix = chart.NormalizedValues; % Get the normalized confusion matrix
	accuracy(ind) = sum(diag(confusionMatrix)) / sum(confusionMatrix(:)); % Calculate accuracy
	title(sprintf('Accuracy = %0.2f%%', accuracy(ind)*100))

	% Set up struct to save data
	neuron_rate_F0(index).putative = nat_data(index).putative;
	neuron_rate_F0(index).CF = nat_data(index).CF;
	neuron_rate_F0(index).MTF = nat_data(index).MTF;
	neuron_rate_F0(index).rate_rep = nat_data(index).oboe_raterep;
	neuron_rate_F0(index).rate = nat_data(index).oboe_rate;
	neuron_rate_F0(index).rate_std = nat_data(index).oboe_rate_std;
	neuron_rate_F0(index).actual = actual2;
	neuron_rate_F0(index).closest = closest2;
	neuron_rate_F0(index).accuracy = accuracy(ind);
	neuron_rate_F0(index).C = C;

end

%% Plot accuracy of each neuron

figure
histogram(accuracy*100,51)
ylabel('# Neurons')
xlabel('Prediction Accuracy (%)')
title('Prediction of F0 for each neuron using average rate')

mean_all = mean(accuracy, 'all');
fprintf('Mean for all = %0.4f\n', mean_all)


%% Save data 

save('Neuron_Rate_F0_Oboe.mat', "neuron_rate_F0")
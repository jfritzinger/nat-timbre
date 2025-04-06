%% model_neuron_rate_F0

%% Load in data

filepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/model_comparisons';
load(fullfile(filepath, 'Data_NT.mat'), 'nat_data')


%% Shape data into model input

% Find all rows with bassoon in them
sesh = [];
for ii = 1:length(nat_data)
	rate = nat_data(ii).bass_rate;
	if ~isempty(rate)
		sesh = [sesh ii];
	end
end
num_data = length(sesh);

% Get all rates for each repetition for bassoon (one example neuron)
for ind = 1:num_data

	index = sesh(ind);
	data = nat_data(index).bass_raterep';
	avg_rate = nat_data(index).bass_rate;
	rate_std = nat_data(index).bass_rate_std;

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
		for ii = 1:40
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
	figure('Position',[136,782,1085,481])
	tiledlayout(1, 2);
	nexttile
	hold on
	bar(avg_rate)
	errorbar(1:40, avg_rate, rate_std/sqrt(20), 'LineStyle','none', 'Color','k')

	% Plot confusion matrix
	actual2 = reshape(actual,[], 20*40);
	closest2 = reshape(closest, [], 20*40);
	nexttile
	C = confusionmat(actual2, closest2);
	confusionchart(C)

	% Calculate accuracy
	chart = confusionchart(actual2,closest2); % Generate confusion chart
	confusionMatrix = chart.NormalizedValues; % Get the normalized confusion matrix
	accuracy(ind) = sum(diag(confusionMatrix)) / sum(confusionMatrix(:)); % Calculate accuracy
	title(sprintf('Accuracy = %0.2f%%', accuracy(ind)*100))

end

%% Plot accuracy of each neuron

figure
histogram(accuracy*100,51)
ylabel('# Neurons')
xlabel('Prediction Accuracy (%)')
title('Prediction of F0 for each neuron using average rate')

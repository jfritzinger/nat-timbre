%% model_neuron_rate_F0
clear

%% Load in data

filepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/model_comparisons';
load(fullfile(filepath, 'Data_NT.mat'), 'nat_data')

%% Shape data into model input

% Find all rows with bassoon in them
sesh = [];
for ii = 1:length(nat_data)
	rate = nat_data(ii).bass_rate;
	rate2 = nat_data(ii).oboe_rate;
	if ~isempty(rate) && ~isempty(rate2)
		sesh = [sesh ii];
	end
end
num_data = length(sesh);
ind_b = 25:40;
ind_o = [1 3:17];
target = 1;

%% Get all rates for each repetition for bassoon (one example neuron)
for ind = 1:num_data

	index = sesh(ind);
	% figure('Position',[136,782,1085,481])
	% tiledlayout(6, 6);

	for target = 1:16
		data_1 = [nat_data(index).bass_raterep(:, ind_b(target)) ...
			nat_data(index).oboe_raterep(:, ind_o(target))];
		idx = (1:20) + 20*(target-1);
		data(idx, :) = data_1;
	end

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

%% Save data

[base, datapath, savepath, ppi] = getPathsNT();
save(fullfile(base, 'model_comparisons', 'Neuron_Rate_Timbre_All.mat'), ...
	"neuron_rate_timbre")

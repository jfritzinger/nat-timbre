%% model_neuron_rate_F0
clear

%% Get list of all timbre stimuli (bassoon)

for iinstr = 2
	if iinstr == 1
		target = 'Oboe';
	else
		target = 'Bassoon';
	end

	if ismac
		fpath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data';
	else
		fpath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\';
	end
	tuning = readtable(fullfile(fpath, 'Tuning.xlsx')); % Load in tuning

	listing = dir(fullfile(fpath, 'waveforms', '*.wav'));
	target_WAV = arrayfun(@(n) contains(listing(n).name, target), 1:numel(listing), 'UniformOutput', false);
	wav_nums =  find(cell2mat(target_WAV));

	d = dir(fullfile(fpath,'waveforms', '*.wav'));
	all_files = sort({d.name});
	nfiles = length(wav_nums);
	wav_npts = zeros(1,nfiles);
	wav_data = cell(1,nfiles);

	for i = 1:nfiles
		files{1,i} = all_files{wav_nums(i)};
	end

	% Sort by frequency of pitch
	index = [];
	note_names = extractBetween(files, 'ff.','.');
	for ii = 1:nfiles % Find index of each note in tuning spreadsheet
		index(ii) = find(strcmp(note_names(ii), tuning.Note));
	end
	pitch_order = tuning.Frequency(index); % Get freqs of each note
	[~, order] = sort(pitch_order); % Sort freqs

	if iinstr == 1
		files_o = files(order);
		note_names_o = note_names(order);
	else
		files_b = files(order);
		note_names_b = note_names(order);
	end
end
bass_pitch = pitch_order(order);

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


		data = [nat_data(index).bass_raterep(:, ind_b(target)) ...
			nat_data(index).oboe_raterep(:, ind_o(target))];
		avg_rate = [nat_data(index).bass_rate(ind_b(target)) nat_data(index).oboe_rate(ind_o(target))];
		rate_std = [nat_data(index).bass_rate_std(ind_b(target)) nat_data(index).oboe_rate_std(ind_o(target))];

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
		actual2 = reshape(actual,[], 20*2);
		closest2 = reshape(closest, [], 20*2);
		C = confusionmat(actual2, closest2);
		chart = confusionchart(actual2,closest2); % Generate confusion chart
		confusionMatrix = chart.NormalizedValues; % Get the normalized confusion matrix
		accuracy(ind, target) = sum(diag(confusionMatrix)) / sum(confusionMatrix(:)); % Calculate accuracy

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
		% % Plot confusion matrix
		% nexttile
		% confusionchart(C)

	end
end

%% Plot accuracy of each neuron
figure
tiledlayout(4, 4)
for ii = 1:16
	nexttile
	histogram(accuracy(:,ii)*100,21)
	mean_F0 = mean(accuracy(:,ii));
	hold on
	xline(mean_F0*100, 'r', 'LineWidth',2)
	ylabel('# Neurons')
	xlabel('Prediction Accuracy (%)')
	title(['Prediction of instrument, F0=' num2str(round(bass_pitch(ind_b(ii))))])
	endn

mean_all = mean(accuracy, 'all');
fprintf('Mean for all = %0.4f\n', mean_all)



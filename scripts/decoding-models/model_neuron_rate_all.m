%% model_neuron_all

%% Load in spreadsheet 

[base, ~, ~, ~] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')

%% Find all rows with bassoon and oboe in them

[sesh, num_data] = getTimbreSessions(nat_data);

%% Set responses 
% 75 * 20 = 1500 responses
tuning = readtable(fullfile(base, 'Tuning.xlsx')); % Load in tuning

% Get bassoon stimulus
target = 'Bassoon';
listing = dir(fullfile(base, 'waveforms', ['*' target '*.wav']));
files = {listing.name};
note_names = extractBetween(files, 'ff.', '.');
[~, index] = ismember(note_names, tuning.Note);
F0s_b = round(tuning.Frequency(index));
[F0s_b, ~] = sort(F0s_b);

% Get bassoon stimulus
target = 'Oboe';
listing = dir(fullfile(base, 'waveforms', ['*' target '*.wav']));
files = {listing.name};
note_names = extractBetween(files, 'ff.', '.');
[~, index] = ismember(note_names, tuning.Note);
F0s_o = round(tuning.Frequency(index));
[F0s_o, ~] = sort(F0s_o);

% Get into 'Response' 
response_b = cell(75,1);
for ii = 1:length(F0s_b)
	response_b{ii} = ['B_' num2str(F0s_b(ii))];
end
for ii = 1:length(F0s_o)
	response_b{ii+length(F0s_b)} = ['O_' num2str(F0s_o(ii))];
end

response = reshape(repmat(response_b, 1, 20)', 1, []);
response = response';

%% Get all rates for each repetition for bassoon (one example neuron)


for ind = 1:num_data

	% Set up rate data 
	index = sesh(ind);
	data = [nat_data(index).bass_raterep nat_data(index).oboe_raterep];

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
		for ii = 1:length(avg_other_rows)
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
	actual2 = reshape(actual,[], 1500);
	closest2 = reshape(closest, [], 1500);
	C = confusionmat(actual2, closest2);
	%chart = confusionchart(actual2,closest2); % Generate confusion chart
	%confusionMatrix = chart.NormalizedValues; % Get the normalized confusion matrix
	accuracy(ind) = sum(diag(C)) / sum(C(:)); % Calculate accuracy

	% Set up struct to save data
	neuron_rate_all(ind).putative = nat_data(index).putative;
	neuron_rate_all(ind).CF = nat_data(index).CF;
	neuron_rate_all(ind).MTF = nat_data(index).MTF;
	neuron_rate_all(ind).rate_rep = data;
	neuron_rate_all(ind).actual = actual2;
	neuron_rate_all(ind).closest = closest2;
	neuron_rate_all(ind).accuracy = accuracy(ind);
	neuron_rate_all(ind).C = C;

	fprintf('%d/%d, %0.2f%% done!\n', ind, num_data, ind/num_data*100)
end

%% Save data 

[base, datapath, savepath, ppi] = getPathsNT();
save(fullfile(base, 'model_comparisons', 'Neuron_Rate_All.mat'), ...
	"neuron_rate_all")


%% Plot outputs 

figure
histogram(accuracy)
hold on
xline(1/75)

%% Run using classification model 




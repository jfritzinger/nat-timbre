%% model_neuron_all

%% Load in spreadsheet 

[base, datapath, savepath, ppi] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')

%% Find all rows with bassoon and oboe in them

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

%% Run model for each neuron

for ind = 1 %:num_data

	% Set up index 
	index = sesh(ind);

	% Create matrix of all data for model
	h_all_bass = [];
	for itarget = 1:length(F0s_b)

		% Arrange bassoon data for SVM
		spikes = nat_data(index).bass_spikerate{itarget}/1000; % ms
		spikereps = nat_data(index).bass_spikerep{itarget};
		min_dis = 0.25;
		edges = 0:min_dis:300;
		t = 0+min_dis/2:min_dis:300-min_dis/2;
		for irep = 1:20
			h_bass(irep, :) = histcounts(spikes(spikereps==irep), edges);
		end
		h_all_bass = [h_all_bass; h_bass];
	end

	h_all_oboe = [];
	for itarget = 1:length(F0s_o)


		spikes = nat_data(index).oboe_spikerate{itarget}/1000; % ms
		spikereps = nat_data(index).oboe_spikerep{itarget};
		for irep = 1:20
			h_oboe(irep, :) = histcounts(spikes(spikereps==irep), edges);
		end
		h_all_oboe = [h_all_oboe; h_oboe];
	end
	h_all = [h_all_bass; h_all_oboe];
	T = array2table(h_all);
	T.response = response;
	predictors = h_all;


	%% Run model



end

%% Plot outputs 






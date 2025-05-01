clear

%% Load in data

filepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/model_comparisons';
load(fullfile(filepath, 'Data_NT2.mat'), 'nat_data')

%% Get correct output of model
target = 'Bassoon';

if ismac
	fpath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data';
else
	fpath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\';
end
tuning = readtable(fullfile(fpath, 'Tuning.xlsx')); % Load in tuning

% Get bassoon stimulus
listing = dir(fullfile(fpath, 'waveforms', ['*' target '*.wav']));
files = {listing.name};
note_names = extractBetween(files, 'ff.', '.');
[~, index] = ismember(note_names, tuning.Note);
F0s1 = tuning.Frequency(index);
[F0s, order] = sort(F0s1);
%F0s = log10(F0s);

%% Get data

h_all2 = [];
for ineuron = 1:3
	if ineuron == 1
		putative = 'R29_TT2_P3_N03'; % Predictions 44, yesterday 62?
	elseif ineuron == 2
		putative = 'R29_TT3_P5_N05'; % Not predicting well, not sure why
	else
		putative = 'R29_TT3_P2_N06'; % Good

	end
	index = find(strcmp({nat_data.putative}, putative));

	h_all = [];
	for itarget = 1:40
		spikes_bass = nat_data(index).bass_spikerate{itarget}/1000; % ms
		spikereps_bass = nat_data(index).bass_spikerep{itarget};

		% Arrange data for SVM
		min_dis = 0.25;
		edges = 0:min_dis:300;
		t = 0+min_dis/2:min_dis:300-min_dis/2;
		for irep = 1:20
			h_bass(irep, :) = histcounts(spikes_bass(spikereps_bass==irep), edges);
		end
		h_all = [h_all; h_bass];
	end

	h_all2 = [h_all2, h_all];

end

% Put data into table
response = reshape(repmat(F0s, 1, 20)', 1, []);
T = array2table(h_all2);
T.response = response';
predictors = h_all2;


%% Call model

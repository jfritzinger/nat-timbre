%% table_of_F0s
clear 


%% Load in data 

target = 'Oboe';
base = getPathsNT;
tuning = readtable(fullfile(base, 'Tuning.xlsx')); % Load in tuning
listing = dir(fullfile(base, 'waveforms', '*.wav'));
target_WAV = arrayfun(@(n) contains(listing(n).name, target), 1:numel(listing), 'UniformOutput', false);
wav_nums =  find(cell2mat(target_WAV));

d = dir(fullfile(base,'waveforms', '*.wav'));
all_files = sort({d.name});
nfiles = length(wav_nums);
wav_npts = zeros(1,nfiles);
wav_data = cell(1,nfiles);
for i = 1:nfiles
	files{1,i} = all_files{wav_nums(i)};
end


%% Sort by frequency of pitch

index = [];
note_names = extractBetween(files, 'ff.','.');
for ii = 1:nfiles % Find index of each note in tuning spreadsheet
	index(ii) = find(strcmp(note_names(ii), tuning.Note));
end
pitch_order = tuning.Frequency(index); % Get freqs of each note



%% Calculate actual F0 

F0s = NaN(nfiles, 1);
peak_harm = NaN(nfiles, 1);
for iii = 1:nfiles

	% Load in file 
	file = files{iii};
	[stim, fs] = audioread(fullfile(base, 'waveforms', file));

	% Welch's power spectral density estimate
	[pxx, f] = pwelch(stim,[],[],[],fs);
	[pks,locs] = findpeaks(pxx); % Find peaks 
	[~, ind] = max(pks);
	peak_harm(iii) = f(locs(ind));

	% Get the first harmonic as the initial F0 guess 
	pks_dB = 10*log10(pks);
	dB_thresh = -70;
	i = find(pks_dB > dB_thresh, 1);
	F0_actual = f(locs(i));
	harmonics(1) = F0_actual;
	current_F0_guess = F0_actual;

	% Iterate through harmonics to get a guess for F0
	max_freq = 10000;
	i = 2;
	this_f = 0;
	while this_f <= max_freq
		next_F_guess = current_F0_guess*i;
		[~,I(i)] = min(abs(next_F_guess - f(locs)));
		if I(i-1)==I(i)
			harmonics(i) = next_F_guess;
			this_f = next_F_guess;
		else
			harmonics(i) = f(locs(I(i)));
			this_f = f(locs(I(i)));
		end
		current_F0_guess = (current_F0_guess*(i-1) + (harmonics(i)-harmonics(i-1)))/i;
		i = i+1;
	end
	F0s(iii) = current_F0_guess;

end

%% Create table 

for ii = 1:nfiles
	Label{ii,1} = ['nt_o_' num2str(ii)];
end
F0_Approx = pitch_order;
F0_Actual = F0s;
Note_Names = note_names';
nt_o = table(Label, F0_Approx, F0_Actual, Note_Names);

writetable(nt_o, 'nt_o.xlsx')
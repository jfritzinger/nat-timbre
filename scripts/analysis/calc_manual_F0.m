%% calc_F0
clear 

%% Load in bassoon data 
target = 'Bassoon';

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

% Sort by frequency of pitch
index = [];
note_names = extractBetween(files, 'ff.','.');
for ii = 1:nfiles % Find index of each note in tuning spreadsheet
	index(ii) = find(strcmp(note_names(ii), tuning.Note));
end
pitch_order = tuning.Frequency(index); % Get freqs of each note
[~, order] = sort(pitch_order); % Sort freqs

files = files(order);
pitch = pitch_order(order);

%% For one stimulus (F0 = 73 Hz), calculate actual F0 

file = files{5};
F0 = pitch(5);
[stim, fs] = audioread(fullfile(base, 'waveforms', file));
dur = length(stim)/fs;
t = 1:dur/fs-1:dur;


%%

[pxx, f] = pwelch(stim,[],[],[],fs);
[pks,locs] = findpeaks(pxx);

pks_dB = 10*log10(pks);
dB_thresh = -70;
i = find(pks_dB > dB_thresh, 1);
F0_actual = f(locs(i));
harmonics(1) = F0_actual;
current_F0_guess = F0_actual;

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
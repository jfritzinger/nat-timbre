%% Calculating ERR and modulation depth of instrument sounds
clear

%% Get list of all timbre stimuli (bassoon)

[base, ~, ~, ~] = getPathsNT();

target = 'Oboe';
%for iinstr = 1:2
	% if iinstr == 1
	% 	target = 'Oboe';
	% else
	% 	target = 'Bassoon';
	% end

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

	pitches = pitch_order(order);
	files_ordered = files(order);
%end

%% Calculate modulation depth 

figure
tiledlayout(6, 6, "TileSpacing","tight", "Padding","tight")
for ii = 1:nfiles

	[stim, fs] = audioread(fullfile(fpath, 'waveforms', files_ordered{ii}));
	F0 = pitches(ii);
	t = linspace(0, length(stim)/fs, length(stim));

	% Extract envelope using Hilbert transform
	env2 = abs(hilbert(stim)); 
	
	% Compute envelope 
	[env, ~] = envelope(stim); %

	% Modulation depth
	Emax = max(env);
	Emin = min(env);
	modDepth(ii) = (Emax - Emin) / (Emax + Emin);

	% Plot envelope
	nexttile
	plot(t, stim, 'k'); hold on;
	plot(t, env, 'r', 'LineWidth', 0.5);
	plot(t, env2, 'b', 'LineWidth', 0.5);
	hold off;
	title(sprintf('Mod Depth = %0.02f', modDepth(ii)))
	xlim([0.1 0.125])

end

%% 

figure
nexttile
scatter(pitches, ERR)
title('ERR')
set(gca, 'xscale', 'log')

nexttile
scatter(pitches, modDepth, 15, 'filled')
title('Modulation Depth')
set(gca, 'xscale', 'log')

%% Save data 

save(fullfile(base, 'ERR_ModDpth_Oboe.mat'), "ERR", "modDepth", "pitches")
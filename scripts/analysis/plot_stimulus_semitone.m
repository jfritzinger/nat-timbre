%% plot_stimulus_semitone
clear

%% Get list of all timbre stimuli (bassoon)

if ismac
	fpath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data';
else
	fpath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\';
end
tuning = readtable(fullfile(fpath, 'Tuning.xlsx')); % Load in tuning

target = 'Bassoon';
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

files = files(order);
pitch = pitch_order(order);

%% Plot one semitone apart

n = 26;

figure('Position',[92,517,849,395])
tiledlayout(2, 1)

for itype = 1:2

	if itype == 1
		ind = [n, n+1];
	else
		ind = [n, n+3];
	end

	target_files = files(ind);
	F0s = pitch(ind);

	nexttile
	colors = {'#1b9e77', '#d95f02', '#7570b3','#e7298a'};
	for j1 = 1:length(target_files)

		hold on
		instrument = target_files{j1};
		[stim, Fs] = audioread(fullfile(fpath, 'waveforms', instrument));
		F0 = F0s(j1);

		% Calculate spectra
		dist = round(F0/4);
		y2 = fft(stim);
		m = abs(y2);
		mdB = 20*log10(m);
		f = (0:length(y2)-1)*Fs/length(y2);
		mdB(mdB<0) = 0;
		f(f>Fs/2) = [];
		mdB = mdB(1:length(f))';

		% Plot
		[pks, locs] = findpeaks(mdB, 'MinPeakDistance', dist);
		freqs = f(locs);
		plot(freqs, pks, '--', 'LineWidth', 1.5, 'Color',...
			colors{j1});
		for ii = 1:length(pks)
			stem(freqs(ii), pks(ii), 'Marker', 'none', 'LineWidth', ...
				2, 'Color', colors{j1});
		end

		% Figure Properties
		fpeaks = F0s(1) * 2.^((0:10));
		xticks(fpeaks)
		xlim([50 10000])
		yticklabels([])
		set(gca, 'XScale', 'log')
		hold on
		box on
		set(gca,'fontsize',14)
		if j1 == length(target_files)
			xlabel('Frequency (Hz)')
		elseif j1 == round(length(target_files)/2)
			ylabel('Mag. (dB SPL)')
			xticklabels([])
		elseif j1 == 1
			xticklabels([])
			title(groupings{idx})
		else
			xticklabels([])
		end
		if itype == 1 && j1 == 1
			title(sprintf('F0 = %0.0f Hz and semitone above, %0.0f Hz', F0, F0s(2)))
		elseif itype == 2 && j1 == 1
			title(sprintf('F0 = %0.0f Hz and whole tone above, %0.0f Hz', F0, F0s(2)))
		end
		ax = gca;
		ax.XGrid='on';
		ax.XMinorGrid='off';
	end
end


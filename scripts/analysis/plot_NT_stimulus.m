%% plot_NT_stimulus.m
%
% Script that plots out the magnitude spectrum for each oboe or bassoon
% stimulus where stimuli that are octaves are plotted together. Also plots
% out the range of F0s used in bassoon and oboe. 
%
% Author: J. Fritzinger
% Created: 2022-09-15; Last revision: 2024-10-14
%
% -------------------------------------------------------------------------
clear 

%% Get list of all timbre stimuli (bassoon)

if ismac
	fpath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data';
else
	fpath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\';
end
tuning = readtable(fullfile(fpath, 'Tuning.xlsx')); % Load in tuning 

target = 'Oboe';
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


%% Are the stimuli spaced out in 1/6 octaves?

% Calculate even array of frequencies
freq_lo = pitch_order(order(1)); 
freq_hi = pitch_order(order(end)); 
fpeaks = freq_lo * 2.^((0:nfiles)/12);

figure
scatter(pitch_order(order), ones(nfiles,1), 'filled');
set(gca, 'xscale', 'log')
xlim([freq_lo-5 freq_hi+50])
xticks([55 110 220 440 880 1760])
xticks(fpeaks)
ax = gca;
ax.XGrid='on';
ax.XMinorGrid='off';

%% Find octave groupings 

groupings = {'A', 'Bb', 'B', 'C', 'Db', 'D', 'Eb', 'E', 'F', 'Gb', 'G', 'Ab'};
note_names = extractBetween(files, 'ff.','.');
for i = 1:length(note_names)
	if length(note_names{i})==2
	note_names{i} = note_names{i}(1);
	else
		note_names{i} = note_names{i}(1:2);
	end
end

for idx = 1:12

	index = strcmp(groupings{idx}, note_names);
	target_files = files(index);
	F0s = pitch_order(index);

	% Plot all octaves groups
	figure('Position',[92,517,849,395])
	tiledlayout(4, 1, 'TileSpacing','none')
	colors = {'#1b9e77', '#d95f02', '#7570b3','#e7298a'};
	for j1 = 1:length(target_files)
		nexttile
		hold on
		instrument = target_files{j1};
		[stim, Fs] = audioread([fpath instrument]);
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
		%plot(f, mdB, 'LineWidth', 2, 'Color', colors{j1});
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
		set(gca,'fontsize',18)
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
		ax = gca;
		ax.XGrid='on';
		ax.XMinorGrid='off';

	end

	savepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/figures/stimuli';
	%print('-dpng', fullfile(savepath, [groupings{idx} '.png']))
	print('-dpng', fullfile(savepath, [num2str(round(F0s(1))) '.png']))

end
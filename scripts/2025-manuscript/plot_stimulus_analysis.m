%% plot_overlapping_stimuli.m
%
% Script that plots the oboe and bassoon sounds that have overlapping F0s. 
%
% Author: J. Fritzinger
% Created: 2022-10-14; Last revision: 2024-10-14
%
% -------------------------------------------------------------------------
clear 

%% Set up figure

addpath('/Users/jfritzinger/Projects/nat-timbre/scripts/helper-functions', '-end')
[base, ~, ~, ppi] = getPathsNT();
figure('Position',[68,278,4.567*ppi,2.7*ppi]);
titlesize = 9;
fontsize = 8;
labelsize = 12;
legsize = 7;
linewidth = 1;

%% Get list of all timbre stimuli (bassoon)

for iinstr = 1:2
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

%% Extract overlapping F0s/stimuli numbers (Bb3 to D5) 
% 16 overlapping 
n_overlap = 16;
lower_lim = 'Bb3';
high_lim = 'D5';

% Bassoon index 
ind_l = find(strcmp(lower_lim, note_names_b));
ind_h = find(strcmp(high_lim, note_names_b));
bassoon_files = files_b(ind_l:ind_h);
F0s = pitch_order(order(ind_l:ind_h));

% Oboe index 
ind_l = find(strcmp(lower_lim, note_names_o));
ind_h = find(strcmp(high_lim, note_names_o));
oboe_files = files_o(ind_l:ind_h);
oboe_files(2) = [];


%% Load and plot all 

colors = {"#0072BD", "#D95319"};
iplot = 1;
for j1 = [2, 4, 9, 11]

		% Load in data 
		[oboe, ~] = audioread(fullfile(fpath, 'waveforms', oboe_files{j1}));
		[bassoon, Fs] = audioread(fullfile(fpath, 'waveforms', bassoon_files{j1}));
		name =  extractBetween(oboe_files{j1}, 'ff.','.');

		% Analyze spectrograms & plot
		h(iplot) = subplot(2, 3, iplot);
		iplot = iplot + 1;
		hold on
		for ii = 1:2

			% Calculate spectra
			F0 = F0s(j1);
			dist = round(F0/4);
			if ii == 1
				stim = oboe;
			else
				stim = bassoon;
			end
			y2 = fft(stim);
			m = abs(y2);
			mdB = 20*log10(m);
			f = (0:length(y2)-1)*Fs/length(y2);
			mdB(mdB<0) = 0;
			f(f>Fs/2) = [];
			mdB = mdB(1:length(f))';

			log_f = log10(f(:,2:end));
			center_of_mass = 10.^((sum(log_f.*mdB(:,2:end), 2))./sum(mdB(:,2:end), 2));

			% Plot
			plot(f/1000, mdB, 'LineWidth',linewidth, 'Color',...
				colors{ii});
			[pks, locs] = findpeaks(mdB, 'MinPeakDistance', dist);
			freqs = f(locs);
			if freqs(1)<F0-10
				plot(freqs(2:end)./1000, pks(2:end), '--', 'LineWidth', linewidth, 'Color',...
					colors{ii});
			else
				plot(freqs./1000, pks, '--', 'LineWidth', linewidth, 'Color',...
					colors{ii});
			end
		end
		
		% Figure Properties	
		xlim([150 8000]./1000)	
		xticks([0.2 0.5 1 2 5 10])
		set(gca, 'XScale', 'log')
		hold on
		box on
		ax = gca;
		ax.XGrid='on';
		ax.XMinorGrid='off';

		if ismember(j1, [9, 11, 13:16])
			xlabel('Frequency (kHz)')
		else
			xticklabels([])
		end

		if ismember(j1, [1, 2, 5, 9, 13])
			ylabel('Mag. (dB SPL)')
		else
			yticklabels([])
		end

		if j1 == 1
			legend('Oboe', '','Bassoon')
		end

		text(0.05, 0.95, name{1}, 'Units', 'normalized', ...
			'VerticalAlignment', 'top', 'FontSize',fontsize)
		set(gca,'fontsize',fontsize)
end

%% 

h(5) = subplot(2, 3, 5);
colors = {"#0072BD", "#D95319", "#EDB120", 	"#7E2F8E", 	"#77AC30", "#4DBEEE"};
targets = {'Oboe', 'Bassoon'};
for iinstr = 1:2
	target = targets{iinstr};

	if ismac
		fpath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data';
	else
		fpath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\';
	end
	listing = dir(fullfile(fpath, 'waveforms', '*.wav'));
	target_WAV = arrayfun(@(n) contains(listing(n).name, target), 1:numel(listing), 'UniformOutput', false);
	wav_nums =  find(cell2mat(target_WAV));
	tuning = readtable(fullfile(fpath, 'Tuning.xlsx'));
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
	for ii = 1:nfiles
		[y(ii,:), Fs] = audioread(fullfile(fpath, 'waveforms', files{ii}));
	end
	t = linspace(0, length(y)/Fs, length(y));
	F0s = pitch_order(order);

	% Analysis 
	for ii = 1:nfiles

		y2 = fft(y(ii, :));
		m = abs(y2);
		mdB = 20*log10(m);
		f = (0:length(y2)-1)*Fs/length(y2);
		mdB(mdB<0) = 0;
		f(f>Fs/2) = [];
		mdB = mdB(1:length(f));
		log_f = log10(f(:,2:end));
		center_of_mass(ii) = 10.^((sum(log_f.*mdB(:,2:end), 2))./sum(mdB(:,2:end), 2));

	end

	% Plot SC vs F0 
	hold on
	scatter(F0s, center_of_mass, 10, 'filled', 'MarkerEdgeColor','k',...
		'MarkerFaceColor',colors{iinstr})
	curve = polyfit(F0s, center_of_mass', 3);
	y1 = polyval(curve,F0s);
	plot(F0s, y1, 'color', colors{iinstr}, 'DisplayName', target, LineWidth=linewidth)
	plot(F0s, F0s, 'Color',[0.7 0.7 0.7], LineWidth=linewidth)
	xlabel('F0s (kHz)')
	ylabel('SC (kHz)')
	set(gca, 'XScale', 'log', 'YScale', 'log')
	labels = [50 100 200 500 1000 2000 5000 10000];
	xticks(labels)
	yticks(labels)
	xticklabels(labels/1000)
	yticklabels(labels/1000)
	xlim([56 1662])
	ylim([300 5000])
	grid on
	clear F0s center_of_mass
	set(gca, 'FontSize', fontsize)
end
legend('Oboe', '', '', 'Bassoon', '', '', 'Location','best', ...
	'fontsize', legsize, 'position', [0.7164,0.2969,0.2082,0.1056])

%% Arrange plots 

left = [0.1 0.37 0.73];
bottom = [0.13 0.555];
width = 0.25;
height = 0.375;

set(h(1), 'position', [left(1) bottom(2) width height])
set(h(2), 'position', [left(2) bottom(2) width height])
set(h(3), 'position', [left(1) bottom(1) width height])
set(h(4), 'position', [left(2) bottom(1) width height])
set(h(5), 'position', [left(3) bottom(2) width height])


left = [0, 0.63];
annotation('textbox',[left(1) 0.98 0.0826 0.0385],'String','A',...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[left(2) 0.98 0.0826 0.0385],'String','B',...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
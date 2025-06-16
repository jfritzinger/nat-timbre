%% plot_overlapping_stimuli.m
%
% Script that plots the oboe and bassoon sounds that have overlapping F0s. 
%
% Author: J. Fritzinger
% Created: 2022-10-14; Last revision: 2024-10-14
%
% -------------------------------------------------------------------------
clear 
save_fig = 1;

%% Set up figure

addpath('/Users/jfritzinger/Projects/nat-timbre/scripts/helper-functions', '-end')
[base, ~, ~, ppi] = getPathsNT();
figure('Position',[68,278,8.6*ppi,5*ppi]);
fontsize = 18;
titlesize = 20;
legsize = 20;
linewidth = 1;
scattersize = 40;

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
		ax = gca;
		ax.XGrid='on';
		ax.XMinorGrid='off';
		yticks([0 25 50])

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
	scatter(F0s, center_of_mass, scattersize, 'filled', 'MarkerEdgeColor','k',...
		'MarkerFaceColor',colors{iinstr})
	% Model 1
	%curve = polyfit(F0s, center_of_mass', 3);
	%y1 = polyval(curve,F0s);
	%plot(F0s, y1, 'color', colors{iinstr}, 'DisplayName', target, LineWidth=linewidth)

	% Model 2
	% mdl = fitlm(F0s, center_of_mass');
	% y1 = mdl.Coefficients{2, 1}*F0s + mdl.Coefficients{1, 1};
	% plot(F0s, y1, 'color', colors{iinstr})

	% Model 3
	mdl = fitlm(log10(F0s), log10(center_of_mass));
	p(1) = mdl.Coefficients.Estimate(2,1);
	p(2) = mdl.Coefficients.Estimate(1,1);
	p(3) = mdl.Coefficients.pValue(2);
	y1 = 10^p(2) .* F0s.^p(1);
	plot(F0s, y1, 'color', colors{iinstr})
	if iinstr == 1
		mdl1msg = sprintf('10^{%0.2f} * x^{%0.2f}, p = %0.4f', p(2), p(1), p(3));
	else 
		mdl2msg = sprintf('10^{%0.2f} * x^{%0.2f}, p = %0.4f', p(2), p(1), p(3));
	end

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
% hleg = legend('Oboe', mdl1msg, '', 'Bassoon', mdl2msg, '', 'Location','best', ...
% 	'fontsize', legsize, 'position', [0.7164,0.2969,0.2082,0.1056]);
hleg = legend('Oboe', '', '', 'Bassoon', '', '', 'Location','best', ...
	'fontsize', legsize, 'box', 'off', 'position', [0.713788425761436,0.329270833333333,0.127221324717286,0.109375]);
hleg.ItemTokenSize = [20, 20];

%% Arrange plots 

left = [0.1 0.37 0.73];
bottom = [0.175 0.59];
width = 0.25;
height = 0.375;

set(h(1), 'position', [left(1) bottom(2) width height])
set(h(2), 'position', [left(2) bottom(2) width height])
set(h(3), 'position', [left(1) bottom(1) width height])
set(h(4), 'position', [left(2) bottom(1) width height])
set(h(5), 'position', [left(3) 0.6 width 0.33])

%% Save figure 

if save_fig == 1
	filename = 'intro_stimulus';
	save_figure_MARC(filename)
end
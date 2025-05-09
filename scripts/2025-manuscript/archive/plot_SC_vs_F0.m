%% plot_SC_vs_F0.m
%
% Script that plots oboe and bassoon stimuli scatter plot of spectral
% centroid vs F0
%
% Author: J. Fritzinger
% Created: 2022-10-14; Last revision: 2024-10-14
%
% -------------------------------------------------------------------------
clear

%% Plot 

figure('Position',[560,507,415,341])
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
	scatter(F0s, center_of_mass, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceColor',colors{iinstr})
	curve = polyfit(F0s, center_of_mass', 3);
	y1 = polyval(curve,F0s);
	plot(F0s, y1, 'color', colors{iinstr}, 'DisplayName', target)
	plot(F0s, F0s, 'Color',[0.7 0.7 0.7])
	xlabel('F0s')
	ylabel('Spectral Centroid (Hz)')
	set(gca, 'XScale', 'log', 'YScale', 'log')
	%legend({'Oboe', '', '', 'Bassoon'}, 'Location','northwest')
	xticks([50 100 200 500 1000 2000 5000 10000])
	yticks([50 100 200 500 1000 2000 5000 10000])
	xlim([56 1662])
	ylim([300 5000])
	grid on
	clear F0s center_of_mass
	set(gca, 'FontSize', 14)
	title('Spectral Centroid vs F0')
end
legend('Oboe', '', '', 'Bassoon', '', '', 'Location','best')

%% Save figure

savepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/figures/manuscript';
exportgraphics(gcf, fullfile(savepath, 'spectral_centroid_v_F0.png'), 'Resolution', 600)


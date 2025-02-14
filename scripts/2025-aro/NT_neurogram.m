%% NT_neurograms 

%% compare_neurograms_stim
clear

%% Sort stimuli

if ismac
	fpath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data';
else
	fpath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\';
end
tuning = readtable(fullfile(fpath, 'Tuning.xlsx')); % Load in tuning
target = 'Bassoon';

% Get all .wav files containing the target instrument name
listing = dir(fullfile(fpath, 'waveforms', ['*' target '*.wav']));
files = {listing.name};

% Extract note names and find corresponding frequencies
note_names = extractBetween(files, 'ff.', '.');
[~, index] = ismember(note_names, tuning.Note);
F0s = tuning.Frequency(index);

% Sort files and frequencies by pitch
[F0s, order] = sort(F0s);
files = files(order);

% Initialize variables for later use (if needed)
nfiles = numel(files);
wav_npts = zeros(1, nfiles);
wav_data = cell(1, nfiles);

if ismac
	savepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/AN_neurogram';
else
	savepath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\AN_neurogram';
end
ii = 25;
F0 = F0s(ii);

%% Plot neurograms
fontsize = 20;

ii = 25;
target = extractBefore(files{ii}, '.');
target_F0 = extractBetween(files{ii}, 'ff.', '.wav');

for iplot = 2
	switch iplot
		case 1
			savefile = sprintf('%s_F0_%s_AN.mat', target, target_F0{1});
			load(fullfile(savepath, savefile),'AN', 'model_params', 'params')
			stim = params.stimulus(ii,:);
		case 2
			savefile = sprintf('%s_F0_%s_AN_TC.mat', target, target_F0{1});
			load(fullfile(savepath, savefile),'AN', 'model_params', 'params')
			stim = params.stim;
		case 3
			savefile = sprintf('%s_F0_%s_AN_ST.mat', target, target_F0{1});
			load(fullfile(savepath, savefile),'AN', 'model_params', 'params')
			stim = params.stim;
	end

	% Analyze AN_sout
	an_sout = squeeze(AN.an_sout);
	CFs = AN.CFs;
	t = linspace(0, params.dur, size(an_sout,2));
	F0 = F0s(ii);
	%NHN = CFs/F0; % Calculate neural harmonic number (CF/F0)
	ylimits = [min(CFs) 5325]; %[min(CFs)/F0 25]; % max(CFs)/F0];
	period_lim = 1/F0*5;

	% Analyze spectrogram
	Fs = params.Fs;
	y2 = fft(stim);
	m = abs(y2);
	mdB = 20*log10(m);
	f = (0:length(y2)-1)*Fs/length(y2);
	mdB(mdB<0) = 0;
	f(f>max(CFs)) = [];
	mdB = mdB(1:length(f));
	%NHN_stim = f/F0;

	% Plot MASD, exclude onset
	t_lims = [0.025 0.225]*Fs;
	an_sout2 = an_sout(:,t_lims(1):t_lims(2)-1);
	t2 = linspace(0, 0.2, size(an_sout2,2));
	spec_diff = diff(an_sout2, 1);
	spec_abs = abs(spec_diff);
	MASD = trapz(t2, spec_abs, 2);
	MASD_no = trapz(t2, spec_diff, 2);

	% Plot stimulus spectrogram
	figure('Position',[515,474,410,400])
	tiledlayout(1, 5, 'TileSpacing','none')
	ax1 = nexttile();
	plot(f, mdB, 'k')
	view(90,90)
	set(gca, 'XDir','reverse')
	ylabel('Neural Harmonic Number (CF/F0)')
	xlim(ylimits)
	ylim([0 70])
	xticks([100 200 500 1000 2000 5000 10000])
	ylabel('Magnitude')
	xlabel('CFs (Hz)')
	set(gca, 'FontSize', fontsize)
	set(gca, 'XScale', 'log')

	% Plot neurogram
	ax2 = nexttile([1,4]);
	s = surf(t2, CFs, an_sout2);
	s.EdgeColor = 'none';
	view(0,90)
	xlim([0 period_lim])
	if iplot == 1
		title('Instrumental Timbre')
	elseif iplot == 2
		title('Zero-Phase Tone Complex')
	else
		title('Zero-Phase Synthetic Timbre')
	end
	xlabel('Time (s)')
	ylim(ylimits)
	set(gca, 'YScale', 'log')
	set(gca, 'FontSize', fontsize)
	yticks([100 200 500 1000 2000 5000 10000])
	yticklabels([])

	%linkaxes([ax1 ax2 ax3],'x')
end

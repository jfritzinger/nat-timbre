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

%% Run models for zero phase flat spectrum complex

% Create stimulus
params.dur = 0.3;
params.Fs = 100000;
params.num_comps = 43;
params.ramp_dur = 0.02;
params.stimdB = 73;

% Get the peak (or center, or edge) freq for this presentation.
npts = params.dur * params.Fs;
Pref = 20e-6; % 20 micropascals
stim = zeros(npts, 1);

t = 0:(1/params.Fs):(params.dur - 1/params.Fs);
comp_freqs = F0 + F0 * (0:(params.num_comps-1)); % frequency comps, centered at this_fpeak, spaced by Delta_F

for icomp = 1:params.num_comps
	fcomp =  comp_freqs(icomp); % frequency of this component
	stim = stim + (20e-6 * 10.^(params.stimdB/20) * sin(2*pi*fcomp*t))'; % scale each compoennt to desired dB SPL
end

% Multiply each col of sin_matrix by appropriate amplitude, add up the
% columns, and multiply by ramp envelope.
stim = Pref * 10.^(params.stimdB/20) * stim / rms(stim); % scale overall RMS level into desired pascals
params.stim = stim .* tukeywin(length(stim),2*params.ramp_dur/params.dur); % apply ramp
params.stim = params.stim';
params.num_stim = 1;

% Plot new stimulus
figure
hold on
dist = round(F0/4);
y2 = fft(params.stim);
m = abs(y2);
mdB = 20*log10(m);
Fs = 100000;
f = (0:length(y2)-1)*Fs/length(y2);
mdB(mdB<0) = 0;
f(f>Fs/2) = [];
mdB = mdB(1:length(f))';
plot(f/1000, mdB, 'LineWidth', 1.5);
xlim([150 8000]./1000)
xticks([0.2 0.5 1 2 5 10])
set(gca, 'XScale', 'log')

% Model parameters
CFs = logspace(log10(125), log10(10000), 300);
model_params.type = 'SFIE';
model_params.range = 1; % 1 = population model, 2 = single cell model
model_params.species = 2; % 1 = cat, 2 = human
model_params.BMF = 100;
model_params.CFs = CFs;
model_params.nAN_fibers_per_CF = 10;
model_params.cohc = 1; % (0-1 where 1 is normal)
model_params.cihc = 1; % (0-1 where 1 is normal)
model_params.nrep = 10; % how many times to run the AN model
model_params.implnt = 1; % 0 = approximate model, 1=exact powerlaw implementation(See Zilany etal., 2009)
model_params.noiseType = 1; % 0 = fixed fGn, 1 = variable fGn) - this is the 'noise' associated with spont. activity of AN fibers - see Zilany et al., 2009. "0" lets you "freeze" it.
model_params.which_IC = 1; % 2 = ModFilt; 1 = SFIE model
model_params.onsetWin = 0.020; % exclusion of onset response, e.g. to omit 1st 50 ms, use 0.050
model_params.fiberType = 3; % AN fiber type. (1 = low SR, 2 = medium SR, 3 = high SR)

% Get single stimulus
ii = 25;
target = extractBefore(files{ii}, '.');
target_F0 = extractBetween(files{ii}, 'ff.', '.wav');
savefile = sprintf('%s_F0_%s_AN_TC.mat', target, target_F0{1});

AN = modelAN(params, model_params); % HSR for IC input
save(fullfile(savepath, savefile),'AN')

elapsedTime = toc(timerVal)/60;
disp(['Model took ' num2str(elapsedTime) ' minutes'])


%% Run model for zero phase complex
clear params 

[stim, Fs] = audioread(fullfile(fpath, 'waveforms', files{ii}));
name = extractBetween(files{ii}, 'ff.','.');

figure
hold on
dist = round(F0/4);
y2 = fft(stim);
m = abs(y2);
mdB = 20*log10(m);
f = (0:length(y2)-1)*Fs/length(y2);
mdB(mdB<0) = 0;
f(f>Fs/2) = [];
mdB = mdB(1:length(f))';

% Plot
plot(f/1000, mdB, 'LineWidth', 1.5);
[pks, locs] = findpeaks(mdB, 'MinPeakDistance', dist);
freqs = f(locs);
plot(freqs./1000, pks, '--', 'LineWidth', 1.5);
xlim([150 8000]./1000)
xticks([0.2 0.5 1 2 5 10])
set(gca, 'XScale', 'log')
hold on
xlabel('Frequency (kHz)')
ylabel('Mag. (dB SPL)')

% Create ST stimulus
params.dur = 0.3;
params.Fs = 100000;
params.spl = 73;
params.ramp_dur = 0.02;

% Time vectors
gate = tukeywin(npts,2*params.ramp_dur/params.dur); %raised cosine ramps
npts = params.dur * params.Fs; % # pts in stimulus
t = (0:(npts-1))/Fs; % time vector
interval = zeros(1,length(t));
harmonics = F0:F0:5350; % component freqs for the central stimulus, when this_fpeak = CF
num_harmonics = length(harmonics);

% Make the stimulus
pks = pks(2:end);
pks_linear = 20e-6 * 10.^(pks/20);
for iharm = 1:num_harmonics
	comp_freq = harmonics(iharm);
	interval = interval + pks_linear(iharm) * sin(2*pi*comp_freq*t);
end
Level_scale = 20e-6*10.^(params.spl/20) * (1/rms(interval)); % overall lienar scalar to bring this centered stimulus up to stimdB
interval = Level_scale * interval; % include dB scaling into the set of harmonic component scalars

stim = interval';
stim = stim.*gate;
params.stim = stim';
params.num_stim = 1;

% Plot new stimulus
hold on
dist = round(F0/4);
y2 = fft(interval);
m = abs(y2);
mdB = 20*log10(m);
f = (0:length(y2)-1)*Fs/length(y2);
mdB(mdB<0) = 0;
f(f>Fs/2) = [];
mdB = mdB(1:length(f))';
plot(f/1000, mdB, 'LineWidth', 1.5);
xlim([150 8000]./1000)
xticks([0.2 0.5 1 2 5 10])
set(gca, 'XScale', 'log')

% Model parameters
CFs = logspace(log10(125), log10(10000), 300);
model_params.type = 'SFIE';
model_params.range = 1; % 1 = population model, 2 = single cell model
model_params.species = 2; % 1 = cat, 2 = human
model_params.BMF = 100;
model_params.CFs = CFs;
model_params.nAN_fibers_per_CF = 10;
model_params.cohc = 1; % (0-1 where 1 is normal)
model_params.cihc = 1; % (0-1 where 1 is normal)
model_params.nrep = 10; % how many times to run the AN model
model_params.implnt = 1; % 0 = approximate model, 1=exact powerlaw implementation(See Zilany etal., 2009)
model_params.noiseType = 1; % 0 = fixed fGn, 1 = variable fGn) - this is the 'noise' associated with spont. activity of AN fibers - see Zilany et al., 2009. "0" lets you "freeze" it.
model_params.which_IC = 1; % 2 = ModFilt; 1 = SFIE model
model_params.onsetWin = 0.020; % exclusion of onset response, e.g. to omit 1st 50 ms, use 0.050
model_params.fiberType = 3; % AN fiber type. (1 = low SR, 2 = medium SR, 3 = high SR)
nstim = 1;

% Get single stimulus
ii = 25;
target = extractBefore(files{ii}, '.');
target_F0 = extractBetween(files{ii}, 'ff.', '.wav');
savefile = sprintf('%s_F0_%s_AN_ST.mat', target, target_F0{1});

AN = modelAN(params, model_params); % HSR for IC input
save(fullfile(savepath, savefile),'AN')

elapsedTime = toc(timerVal)/60;
disp(['Model took ' num2str(elapsedTime) ' minutes'])


%% Plot neurograms

for iplot = 1:3
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
	ylimits = [min(CFs)/F0 25]; % max(CFs)/F0];
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
	NHN_stim = f/F0;

	% Plot MASD, exclude onset
	t_lims = [0.025 0.225]*Fs;
	an_sout2 = an_sout(:,t_lims(1):t_lims(2)-1);
	t2 = linspace(0, 0.2, size(an_sout2,2));
	spec_diff = diff(an_sout2, 1);
	%spec_abs = abs(spec_diff);
	MASD = trapz(t2, spec_diff, 2);

	% Plot stimulus spectrogram
	figure('Position',[515,275,1082,573])
	tiledlayout(1, 7, 'TileSpacing','none')
	ax1 = nexttile([1,3]);
	plot(NHN_stim, mdB, 'k')
	view(90,90)
	set(gca, 'XDir','reverse')
	title('Stimulus Spectrogram')
	ylabel('Neural Harmonic Number (CF/F0)')
	xlim(ylimits)
	ylabel('Magnitude')
	set(gca, 'FontSize', 16)

	% Plot neurogram
	ax2 = nexttile([1,3]);
	imagesc(NHN, t, an_sout2')
	set(gca, 'YDir', 'normal');
	set(gca, 'XDir', 'reverse');
	view(90,90)
	ylim([0 period_lim])
	title('AN Neurogram')
	ylabel('Time (s)')
	xlim(ylimits)
	set(gca, 'FontSize', 16)

	% Plot MASD
	ax3 = nexttile;
	plot(NHN(1:end-1), MASD, 'k')
	set(gca, 'XDir','reverse')
	view(90,90)
	hold on
	xticks([])
	title('MASD')
	set(gca,'FontSize',16)
	xlim(ylimits)

	linkaxes([ax1 ax2 ax3],'x')
end

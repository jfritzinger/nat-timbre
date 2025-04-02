%% plot_model_temporal 
clear 

%% Parameters
CF = 1200; 

% Stimulus parameters
params.Fs = 100000;
params.nrep = 1;
params.target = 'Bassoon';
params.signal_onset_delay = 0;
params.noise_ramp_dur = 0.002;
params.signal_spls = 73;
params.reptim = 0.6;
params.noise_state = 0;
params.noise_shape = [];
params.mnrep = 1;
params.noise_spls = 0;
params.signal_ramp_dur = 0.002;
params = generate_NT(params);
params.num_stim = length(params.files);

% Model parameters
model_params.type = 'SFIE';
model_params.range = 2; % 1 = population model, 2 = single cell model
model_params.species = 1; % 1 = cat, 2 = human
model_params.BMF = 100;
model_params.CF_range = 1200;
model_params.num_CFs = 1;
model_params.CFs = 1200;
model_params.nAN_fibers_per_CF = 5;
model_params.cohc = 1; % (0-1 where 1 is normal)
model_params.cihc = 1; % (0-1 where 1 is normal)
model_params.nrep = 10; % how many times to run the AN model
model_params.implnt = 1; % 0 = approximate model, 1=exact powerlaw 
% implementation(See Zilany etal., 2009)
model_params.noiseType = 1; % 0 = fixed fGn, 1 = variable fGn) - 
% this is the 'noise' associated with spont. activity of AN fibers - 
% see Zilany et al., 2009. "0" lets you "freeze" it.
model_params.which_IC = 1; % 2 = ModFilt; 1 = SFIE model
model_params.onsetWin = 0.020; % exclusion of onset response, e.g. to 
% omit 1st 50 ms, use 0.050
model_params.fiberType = 3; % AN fiber type. (1 = low SR, 2 = medium 
% SR, 3 = high SR)
model_params.Fs = 100000;


%% Model 

AN_HSR = modelAN(params, model_params); % HSR for IC input
SFIE = wrapperIC(AN_HSR.an_sout, params, model_params); % SFIE output
% save('/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Aim 2 - 
% Timbre/Data/Intro_ModelResponse.mat', 'AN_HSR', 'SFIE')

% addpath '/Users/jfritzinger/Projects/synth-timbre/scripts/helper-functions'
% [base, datapath, ~] = getPaths();
% filename = 'Intro_ModelResponse.mat';
% load(fullfile(datapath, filename), 'AN_HSR', 'SFIE')
avAN = AN_HSR.average_AN_sout;
avBE = SFIE.average_ic_sout_BE;
avBS = SFIE.average_ic_sout_BS;


%% AN Plot
figure
tiledlayout(3,1)
font_size = 12;
pitches = sort(params.pitch_order);

nexttile
hold on
[rate, ~, pitch] = plotNT(params, avAN, 0);
plot(pitches, rate, 'linewidth', 2, 'color', '#117733');
set(gca,'fontsize',font_size)
grid on
set(gca, 'XScale', 'log')
ylim([0 300])
xlabel('F0 (Hz)')
ylabel('Avg. Rate (sp/s)')
title('AN Avg. Rate')
xlim([55 600])
xticks([55 110 220 440 600])

%% IC BE Plot
nexttile
yyaxis left
hold on
[rate, ~, pitch] = plotNT(params, avBE, 0);
plot(pitches, rate, 'linewidth', 2, 'color', [0, 0.4470, 0.7410]);
set(gca,'fontsize',font_size)
grid on
set(gca, 'XScale', 'log')
ylim([0 50])
yticks([0 30 60])
title('BE Avg. Rate')
ylabel('Avg. Rate (sp/s)')
xlabel('F0 (Hz)')
xticks([55 110 220 440 600])

% IC temporal plot 
fs = params.Fs;
spike_hist = squeeze(SFIE.ic_BE);
VS = calcVS(params, spike_hist, fs, pitches);
yyaxis right
plot(pitches, VS, 'linewidth', 2, 'color', '#d95f02')
ylabel('|fft| @ 200 Hz')
set(gca,'fontsize',font_size)
title('BE: fft @ 200 Hz')
xlabel('F0 (Hz)')
xlim([55 600])
xticks([55 110 220 440 600])

%% IC BS Plot
nexttile
yyaxis left
hold on
[rate, ~, pitch] = plotNT(params, avBS, 0);
plot(pitches, rate, 'linewidth', 2, 'color', [0, 0.4470, 0.7410]);
set(gca,'fontsize',font_size)
grid on
set(gca, 'XScale', 'log')
ylim([0 50])
yticks([0 30 60])
title('BE Avg. Rate')
ylabel('Avg. Rate (sp/s)')
xlabel('pitches (Hz)')
xticks([55 110 220 440 600])

% IC temporal plot 
fs = params.Fs;
spike_hist = squeeze(SFIE.ic_BS);
VS = calcVS(params, spike_hist, fs, pitches);
yyaxis right
plot(pitches, VS, 'linewidth', 2, 'color', '#d95f02')
ylabel('|fft| @ 200 Hz')
set(gca,'fontsize',font_size)
title('BE: fft @ 200 Hz')
xlabel('F0 (Hz)')
xlim([55 600])
xticks([55 110 220 440 600])

%% Calculate PSTHs and period histograms for model PSTH 

figure('Position',[1700,173,463,1048])
tiledlayout(1, 6, 'TileSpacing','compact', 'Padding','compact')
num_stim = length(pitches);
for ii = 1:2
	if ii == 1
		spike_hist = squeeze(SFIE.ic_BE);
		name = 'BE';
	else
		spike_hist = squeeze(SFIE.ic_BS);
		name = 'BS';
	end

	nexttile([1, 2])
	hold on
	max_rate = max(spike_hist, [], 'all')/2;
	onsetwin = 0.05;
	for j = 1:num_stim

		% Plot PSTHs
		spike_wo_onset = spike_hist(j, onsetwin*fs:params.dur*fs-1);
		t = linspace(0.05, 0.3, fs*0.25);
		offset = max_rate * (j-1);
		plot(t, spike_wo_onset + offset);
	end
	ylim([0 max_rate*num_stim])
	xlabel('Time (s)')
	box on
	yticks(linspace(max_rate/2, max_rate*100-max_rate/2, 100))
	yticklabels(pitches)
	xlim([0.05 0.3])
	grid on
	title([name ', PSTH (excluding onset)'])
	set(gca, 'fontsize',12)

	% Calculate period histogram
	nexttile
	max_rate = max(spike_hist, [], 'all')/2;
	onsetwin = 0.05;
	hold on
	for j = 1:num_stim

		% Plot PSTHs
		
		freq = pitches(j); % Stimulus frequency in Hz
		period = 1 / freq; % Period in ms
		samples_per_period = floor(fs*period);

		num_periods = floor(0.25/period);
		t = linspace(0.05, 0.05+num_periods*period, fs*num_periods*period);
		spike_wo_onset = spike_hist(j, onsetwin*fs:floor(0.05*fs+samples_per_period*num_periods)-1);

		period_hist = reshape(spike_wo_onset, samples_per_period,[]);
		avg = mean(period_hist, 2);
		t_period = linspace(0, period, fs*period);

		offset = max_rate * (j-1);
		patch([t_period flip(t_period)],[avg+offset; ...
			repmat(offset,samples_per_period, 1)], 'k')

	end
	ylim([0 max_rate*num_stim])
	xlabel('Time (s)')
	box on
	yticks(linspace(max_rate/2, max_rate*100-max_rate/2, 100))
	yticklabels(pitches)
	xlim([0 0.0175])
	%xticks(0:0.001:0.005)
	grid on
	title([name ', Period Histogram'])
	set(gca, 'fontsize',12)
end

%% FUNCTIONS --------------------------------------------------------------

function R = calcVS(params, spike_hist, fs, pitches)
t = linspace(0, 0.25, fs*0.25);
onsetwin = 0.05; % ms
num_stim = length(pitches);
for ii = 1:num_stim
	f = pitches(ii);
	r = spike_hist(ii, onsetwin*fs:params.dur*fs-1);
	R(ii) = abs(1/sum(r) * sum(r .* exp(1i * 2*pi * f .* t)));
end
end

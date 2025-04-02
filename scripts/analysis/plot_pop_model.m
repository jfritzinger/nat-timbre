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

% Cut to one stimulus, F0 = 
index = find(params.order==find(params.pitch_order==110));
params.stimulus = params.stim;
params.stim = params.stim(index,:);
params.num_stim = 1;

% Model parameters
model_params.type = 'SFIE';
model_params.range = 2; % 1 = population model, 2 = single cell model
model_params.species = 1; % 1 = cat, 2 = human
model_params.BMF = 100;
model_params.CF_range = [125 10000];
model_params.num_CFs = 100;
model_params.CFs = logspace(log10(125), log10(10000), 100);
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
%save('/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Aim 2 - 
% Timbre/Data/Intro_ModelResponse.mat', 'AN_HSR', 'SFIE')

% addpath '/Users/jfritzinger/Projects/synth-timbre/scripts/helper-functions'
% [base, datapath, ~] = getPaths();
% filename = 'Intro_ModelResponse.mat';
% load(fullfile(datapath, filename), 'AN_HSR', 'SFIE')
avAN = AN_HSR.average_AN_sout;
avBE = SFIE.average_ic_sout_BE;
avBS = SFIE.average_ic_sout_BS;
CFs = AN_HSR.CFs;

%% Plot stimulus 
figure('Position',[560,64,560,835])
tiledlayout(6, 1)
font_size = 12;

F0 = params.pitch_order(params.order(index));
dist = round(F0/4);
stim = params.stim;
Fs = params.Fs;
y2 = fft(stim);
m = abs(y2);
mdB = 20*log10(m);
f = (0:length(y2)-1)*Fs/length(y2);
mdB(mdB<0) = 0;
f(f>Fs/2) = [];
mdB = mdB(1:length(f))';
log_f = log10(f(:,2:end));

% Plot
nexttile
hold on
plot(f/1000, mdB, 'LineWidth', 1.5);
[pks, locs] = findpeaks(mdB, 'MinPeakDistance', dist);
freqs = f(locs);
if freqs(1)<F0-10
	plot(freqs(2:end)./1000, pks(2:end), '--', 'LineWidth', 1.5);
else
	plot(freqs./1000, pks, '--', 'LineWidth', 1.5);
end

% Plot envelope of stimulus
plot_range = [0.125 10];
ticks = [0.05 0.1 0.2 0.5 1 2 5 10];
xlim(plot_range)
xticks(ticks)
set(gca,'fontsize',font_size)
ylabel('Level (dB SPL)')
grid on
set(gca, 'XScale', 'log')
xlabel('CF (Hz)')
title(['Stimulus, F0 = ' num2str(F0)])

%% AN Plot
nexttile
hold on
plot(CFs/1000, avAN, 'linewidth', 2, 'color', '#117733');
set(gca,'fontsize',font_size)
grid on
set(gca, 'XScale', 'log')
xlim(plot_range)
ylim([0 300])
xticks(ticks)
xlabel('CF (Hz)')
ylabel('Avg. Rate (sp/s)')
title('AN Avg. Rate')

%% IC BE Plot
nexttile
hold on
plot(CFs/1000, avBE, 'linewidth', 2, 'color', [0, 0.4470, 0.7410]);
set(gca,'fontsize',font_size)
grid on
set(gca, 'XScale', 'log')
xlim(plot_range)
ylim([0 70])
xticks(ticks)
yticks([0 30 60])
xticklabels([])
title('BE Avg. Rate')
ylabel('Avg. Rate (sp/s)')
xlabel('CF (Hz)')

% IC temporal plot 
fs = params.Fs;
spike_hist = squeeze(SFIE.ic_BE);
%VS = calcFFT(spike_hist, params, fs);
VS = calcVS(params, spike_hist, fs, F0);
nexttile
plot(CFs/1000, VS, 'linewidth', 2, 'color', '#d95f02')
xlim(plot_range)
set(gca, 'XScale', 'log')
xticks(ticks)
ylabel(['|fft| @' num2str(F0) ' Hz'])
%ylim([0 100])
set(gca,'fontsize',font_size)
title(['BE: |fft| @' num2str(F0) ' Hz'])
xlabel('CF (Hz)')

%% IC BS Plot
nexttile
hold on
plot(CFs/1000, avBS, 'linewidth', 2, 'color', [0, 0.4470, 0.7410]);
set(gca,'fontsize',font_size)
grid on
set(gca, 'XScale', 'log')
xlim(plot_range)
ylim([0 50])
xticks(ticks)
yticks([0 30 60])
xlabel('CF (Hz)')
title('BS Avg. Rate')
ylabel('Avg. Rate (sp/s)')

% IC temporal plot 
fs = params.Fs;
spike_hist = squeeze(SFIE.ic_BS);
%VS = calcFFT(spike_hist, params, fs);
VS = calcVS(params, spike_hist, fs, F0);
nexttile
plot(CFs/1000, VS, 'linewidth', 2, 'Color','#d95f02')
xlim(plot_range)
set(gca, 'XScale', 'log')
xticks(ticks)
ylabel(['|fft| @' num2str(F0) ' Hz'])
%ylim([0 100])
set(gca,'fontsize',font_size)
title(['BS: |fft| @' num2str(F0) ' Hz'])
xlabel('CF (Hz)')

%% Calculate PSTHs and period histograms for model PSTH 

figure('Position',[1700,173,463,1048])
tiledlayout(1, 6, 'TileSpacing','compact', 'Padding','compact')

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
	for j = 1:100

		% Plot PSTHs
		spike_wo_onset = spike_hist(j, onsetwin*fs:params.dur*fs-1);
		t = linspace(0.05, 0.3, fs*0.25);
		offset = max_rate * (j-1);
		plot(t, spike_wo_onset + offset);
	end
	ylim([0 max_rate*100])
	xlabel('Time (s)')
	box on
	yticks(linspace(max_rate/2, max_rate*100-max_rate/2, 100))
	yticklabels(round(CFs))
	xlim([0.05 0.3])
	grid on
	title([name ', PSTH (excluding onset)'])
	set(gca, 'fontsize',12)

	% Calculate period histogram
	nexttile
	max_rate = max(spike_hist, [], 'all')/2;
	onsetwin = 0.05;
	hold on
	for j = 1:100

		% Plot PSTHs
		freq = F0; % Stimulus frequency in Hz
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
	ylim([0 max_rate*100])
	xlabel('Time (s)')
	box on
	yticks(linspace(max_rate/2, max_rate*100-max_rate/2, 100))
	yticklabels(round(CFs))
	xlim([0 0.005])
	xticks(0:0.001:0.005)
	grid on
	title([name ', Period Histogram'])
	set(gca, 'fontsize',12)
end

%% FUNCTIONS --------------------------------------------------------------

function R = calcVS(params, spike_hist, fs, F0)
t = linspace(0, 0.25, fs*0.25);
f = F0;
onsetwin = 0.05; % ms
for ii = 1:100
	r = spike_hist(ii, onsetwin*fs:params.dur*fs-1);
	R(ii) = abs(1/sum(r) * sum(r .* exp(1i * 2*pi * f .* t)));
end
end

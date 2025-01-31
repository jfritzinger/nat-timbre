%% stimulus_waveforms
clear 

if ismac
	fpath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data';
else
	fpath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\';
end

if ismac
	stim_dir = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/waveforms';
else
	stim_dir = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\waveforms';
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


% Load in tuning sheet
tuning_dir = stim_dir;
tuning = readtable(fullfile(tuning_dir, 'Tuning.xlsx'));


%% Generate stimulus

% % Set up stimuli 
% params.Fs = 100000;
% params.mnrep = 1;
% params.dur = 0.2;
% params.target = 'Bassoon';
% params.signal_onset_delay = 0;
% params.noise_ramp_dur = 0.02;
% params.signal_spls = 73;
% params.nrep = 1;
% params.reptim = 0.6;
% params.noise_state = 0;
% params.noise_shape = [];
% params.noise_spls = 0;
% params.signal_ramp_dur = 0.02;
% [params] = generate_NT(params);
% params.num_stim = size(params.stim, 1);


%% Plot temporal waveforms 
fontsize = 20;

% Load in wav file
target = 'Bassoon.ff.Bb3.wav';
[data,wav_fs] = audioread(fullfile(stim_dir, target));
t = linspace(0, length(data)/wav_fs, length(data));

figure('Position',[496,414,716,425])
tiledlayout(2, 1)
nexttile
plot(t.*1000, data)
set(gca, 'fontsize', fontsize)
yticklabels([])
xlim([0 150])

target = 'Oboe.ff.Bb3.wav';
[data,wav_fs] = audioread(fullfile(stim_dir, target));
nexttile
t = linspace(0, length(data)/wav_fs, length(data));
plot(t.*1000, data)
xlabel('Time (ms)')
set(gca, 'fontsize', fontsize)
yticklabels([])
xlim([0 150])



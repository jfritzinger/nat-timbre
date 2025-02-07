%% stimulus_waveforms
clear 

if ismac
	stim_dir = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/waveforms';
else
	stim_dir = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\waveforms';
end

%% Plot temporal waveforms 
fontsize = 20;

% Load in wav file
target = 'Bassoon.ff.Bb3.wav';
[data,wav_fs] = audioread(fullfile(stim_dir, target));
t = linspace(0, length(data)/wav_fs, length(data));

figure('Position',[144,574,464,265])
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



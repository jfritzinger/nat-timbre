%% stimulus_annotated
clear

if ismac
	stim_dir = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/waveforms';
else
	stim_dir = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\waveforms';
end

%% Plot temporal waveforms
fontsize = 22;

% Load in wav file
target = 'Bassoon.ff.Bb3.wav';
[stim,Fs] = audioread(fullfile(stim_dir, target));
t = linspace(0, length(stim)/Fs, length(stim));
F0 = 233;

nfft = 2^nextpow2(length(stim)*2); % Double the length for better resolution

% Apply Hanning window to reduce spectral leakage
win = hanning(length(stim));
stim_windowed = stim .* win;

y2 = fft(stim_windowed, nfft); % Compute FFT with zero-padding
m = abs(y2); % Calculate magnitude spectrum
mdB = 20*log10(m); % Convert to dB with improved dynamic range
mdB = mdB - max(mdB); % Normalize to 0 dB maximum
f = (0:length(y2)-1)*Fs/length(y2); % Create frequency axis
noise_floor = -60; % Apply noise floor threshold
mdB(mdB < noise_floor) = noise_floor;

% Keep only positive frequencies up to Nyquist
f(f>Fs/2) = [];
mdB = mdB(1:length(f))';

figure1 = figure('Position', [260,456,600,300]);
axes('Position',[0.1 0.2 0.85 0.6])
hold on
plot(f/1000, mdB, 'LineWidth', 3);


dist = round(F0/2);
[pks, locs] = findpeaks(mdB, 'MinPeakDistance', dist); % Adjust prominence threshold
freqs = f(locs);
if freqs(1) < F0-10
	harmonic_freqs = freqs(2:end);
	harmonic_pks = pks(2:end);
else
	harmonic_freqs = freqs;
	harmonic_pks = pks;
end
plot(harmonic_freqs./1000, harmonic_pks, '--', ...
	'MarkerSize', 8, 'LineWidth', 3, 'Color','#31a354');

xlim([150 8000]./1000)
ylim([noise_floor 5])
xticks([0.2 0.5 1 2 5 8])
yticks([-60 -40 -20 0])
yticklabels([0 20 40 60])
set(gca, 'XScale', 'log')
grid on
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.GridAlpha = 0.3;
ax.MinorGridAlpha = 0.2;
xlabel('Frequency (kHz)', 'FontSize', fontsize)
ylabel('Magnitude (dB)', 'FontSize', fontsize)
set(gca, 'FontSize', fontsize)
box on

log_f = log10(f(2:end));
mdB = mdB+60;
center_of_mass = 10.^((sum(log_f.*mdB(2:end)))./sum(mdB(2:end)));
xline(center_of_mass/1000, ':', 'linewidth', 4)

%%

% Create textbox
annotation(figure1,'textbox',...
	[0.0150000000000019 0.872333333333339 0.465833333333335 0.11],...
	'String',{'Fundamental frequency (F0)'},...
	'FontSize',20,...
	'FitBoxToText','off',...
	'EdgeColor','none');
annotation(figure1,'arrow',[0.113333333333334 0.181666666666667],...
	[0.889000000000003 0.610000000000002]);

% Create textbox
annotation(figure1,'textbox',...
	[0.800000000000014 0.535666666666676 0.265833333333334 0.110000000000001],...
	'String',['Spectral',sprintf('\n'),'peaks'],...
	'FontSize',20,...
	'FitBoxToText','off',...
	'EdgeColor','none');
annotation(figure1,'arrow',[0.796666666666667 0.746666666666667],...
	[0.53 0.403333333333333]);
annotation(figure1,'arrow',[0.79 0.658333333333333],...
	[0.582333333333333 0.49]);
annotation(figure1,'arrow',[0.82 0.813333333333333],...
	[0.455666666666667 0.353333333333333]);


% Create textbox
annotation(figure1,'textbox',...
	[0.8 0.749000000000002 0.205 0.11],'String',{'Harmonics'},...
	'FontSize',20,...
	'EdgeColor','none');
annotation(figure1,'line',[0.3 0.928],[0.774 0.774]);
annotation(figure1,'line',[0.3 0.3],[0.682 0.773]);
annotation(figure1,'line',[0.928 0.928],[0.682 0.773]);

% Create textbox
annotation(figure1,'textbox',...
	[0.48 0.792333333333339 0.176666666666667 0.186666666666667],...
	'String',{'Spectral','centroid'},...
	'FontSize',20,...
	'EdgeColor','none');


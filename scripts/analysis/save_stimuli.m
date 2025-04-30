% Script to read in all of the instrument sounds and manually cut them to
% the right length for the natural timbre stimuli
% Cut sounds to 300ms, in the most steady state area
% J. Fritzinger, 7/7/2021

%% Load in names of music files
clear all
close all 

filepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/zipped_waveforms/Bassoon.ff.stereo';
files = dir(fullfile(filepath, '*.aif'));
fpath = files(1).folder;
for ii = 1:size(files, 1)
    filenames{ii} = fullfile(fpath, files(ii).name);
end

% Load in music file 1 at a time 
index = 10;
[y, Fs] = audioread(filenames{index});
t = linspace(0, length(y)/Fs, length(y));

figure;
set(gcf,'color','w');
plot(t, y(:,1))
title('Full Waveform')
soundsc(y(:,1), Fs);

%% Cut and save 

starttime = 0.61;
dist = 240;

dur = 0.3; % ms
t_cut = linspace(0, dur, Fs*dur);
y_cut = y(starttime*Fs:(starttime+dur)*Fs-1, 1);
%soundsc(y_cut, Fs);
filename_new = erase(files(index).name,".stereo.aif");
audiowrite(['Cut Waveforms\' filename_new '.wav'], y_cut, Fs); 

% Save figure associated with the music file
fig = figure();
fig.Position = [415,407,1091,439];
set(gcf,'color','w');
subplot(2, 2, 1)
plot(t, y(:,1))
title('Full Waveform')
xlabel('Time (s)')
ylabel('Amplitude')
subplot(2, 2, 2)
plot(t_cut, y_cut);
title('Cut Waveform')
xlabel('Time (s)')
ylabel('Amplitude')

% % Scale waveforms to same level 
% spl = 70;
% ramp_dur = 0.02;
% gate = tukeywin(length(y_cut), 2*ramp_dur/dur);
% rms_wavfile = rms(y_cut, 2);                       % Current rms
% rms_Pa = (10^(spl/20))*20e-6;                      % Desired rms amplitude
% y_scaled = y_cut.*(rms_Pa./rms_wavfile);           % Signal in Pa
% y_scaled = y_scaled.*gate;                         % Gated signal with raised cosine ramps
% 
% % Plot scaled waveforms 
% subplot(2, 3, 3)
% plot(t_cut, y_scaled)
% title('Scaled and Ramped Waveform')
% xlabel('Time (s)')
% ylabel('Amplitude')

% Calculate spectra
y2 = fft(y_cut);
m = abs(y2);
mdB = 20*log10(m);
f = (0:length(y2)-1)*Fs/length(y2);
mdB(mdB<0) = 0;
f(f>Fs/2) = [];
mdB = mdB(1:length(f))';
log_f = log10(f(:,2:end));
center_of_mass = 10.^((sum(log_f.*mdB(:,2:end), 2))./sum(mdB(:,2:end), 2));

% Plot spectra 
subplot(2, 2, 3:4)
hold on
semilogx(f,mdB);
[pks, locs] = findpeaks(mdB, 'MinPeakDistance', dist);
plot(f(locs),pks,'linewidth', 1.2, 'Color', [150 150 150]/255);
xline(center_of_mass, '--', 'linewidth', 1.5);
xlim([50 14000])
ylim([0 70])
grid on
set(gca, 'XScale', 'log')
legend(['SC = ' num2str(round(center_of_mass)) 'Hz']);
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB SPL)')
xticks([100 200 500 1000 2000 5000 10000 20000])
box on
title(['Spectrum of ' files(index).name])

% Save figure 
%exportgraphics(fig, [files(index).name '.jpg'], 'Resolution', 300);
saveas(fig, ['Figures\' filename_new '.jpg']);

%% NT_spectrograms.m
%
% J. Fritzinger

clear

%% Get list of all timbre stimuli (bassoon)

if ismac
	fpath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data';
else
	fpath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\';
end
tuning = readtable(fullfile(fpath, 'Tuning.xlsx')); % Load in tuning 

target = 'Bassoon';
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


%% 

for ii = 1 %:40

	target_files = files{ii};
	[stim, Fs] = audioread([fpath target_files]);

	% Calculate the spectrogram of a stimulus
	windowLength = 1000; %hanning(256); %50; % Length of the analysis window (in samples)
	noverlap = []; %[]; %49; % Overlap size between consecutive windows (in samples)
	nfft = 512.; % Length of the FFT

	figure('Position',[59,831,861,362])
	tiledlayout(1, 3)
	%nexttile
	[s,f,t] = spectrogram(stim, windowLength, noverlap, nfft, Fs);
	%spectrogram(stim, windowLength, noverlap, nfft, Fs)

	nexttile
	p = abs(s).^2;
	imagesc(t, f, 10*log10(p));
	xlabel('Time')
	ylabel('Frequency (Hz)')
	set(gca,'Ydir','normal')
	maxdB = max(10*log10(p(:)));
	mindB = max(10*log10(p(:))) - 60;  % Set minimum to 60 dB below maximum
	clim([mindB maxdB]);
	ylim([0.1, 10000])
	% %colormap('bone')

	nexttile
	imagesc(t, f, abs(s));
	xlabel('Time')
	ylabel('Frequency (Hz)')
	set(gca,'Ydir','normal')
	ylim([0.1, 10000])
	%clim([-0.6 2.2])
	colorbar

	nexttile
	imagesc(t, f, angle(s));
	xlabel('Time')
	ylabel('Frequency (Hz)')
	set(gca,'Ydir','normal')
	ylim([0.1, 10000])
	%clim([-0.6 2.2])
	colorbar
			

end
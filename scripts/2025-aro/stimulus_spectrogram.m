%% chirp_analysis_wrapper.m
%
% Loads in natural timbre stimuli and runs through chirp_analysis 
%
% 
% 
clear 

%% Get list of all timbre stimuli (bassoon)

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

%% 

ind = 24; % for bassoon 220 Hz 
%ind = 1;
target = extractBefore(files{ind}, '.');
target_F0 = F0s(ind);
file = fullfile(fpath,'waveforms', files{ind});

%% Analysis

% Reads in file
[data,fs] = audioread(file);

% Creates filename for analyzed data
if ismac
	savepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/decomped';
else
	savepath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\decomped';
end
savefile = sprintf('%s_F0_%0.0f_decomp.mat', target, target_F0);
F0 = target_F0;

% Sets time vector & max frequency
dur = length(data)/fs; %in sec
t = (1/fs):(1/fs):dur;
t = t(1:length(data)); % not sure why t and data are different lengths, this remedies that

% Determines maximum frequency to analyze
max_freq = F0*25; % JBF change, originally 3500; 
if max_freq>20000
	max_freq = 20000;
end

temp = load(fullfile(savepath, savefile));
decomp_info = temp.decomp_info;

% get decomp variables out of struct
starti = decomp_info(1).start_inds{1};
endi = decomp_info(1).end_inds{1};
num_parts = decomp_info(1).num_parts;

% construct "not-spectrogram" spectrogram using phase first, append
% waveforms of cosines with appropriate phase into the same array
phase_spect = zeros(length(decomp_info.freq{1}),length(t));
for icomp = 1:length(decomp_info.freq{1})
	A2 = decomp_info.A{1}(icomp);
	f2 = decomp_info.freq{1}(icomp);
	phi2 = decomp_info.phi{1}(icomp);
	phase_spect(icomp,:) = (A2*cos(2*pi*t*f2 + phi2));
end

% manually unwrap phase so that it is monotonically decreasing
for parti = 1:num_parts
	if parti == 1
		data_part = data(starti(parti):endi(parti));
		A = decomp_info(parti).A{:};
		freq = decomp_info(parti).freq{:};
		phi = decomp_info(parti).phi{:};
		uwphi = [];

		pwelch(data_part,[],[],[],fs);
		xlim([0 5])
		for iphi = 1:length(phi)
			if iphi == 1
				uwphi(iphi) = phi(iphi);
			else
				temp = phi(iphi);
				while temp > uwphi(iphi-1)
					temp = temp - 2*pi;
				end
				uwphi(iphi) = temp;
			end
		end
	end

	% Reconstruct signal using these phases
	temp_pin = zeros(1,length(t(starti(parti):endi(parti))));
	for iharm = 1:length(freq)
		% without magnitude
		temp_pin = temp_pin + cos(2*pi*t(starti(parti):endi(parti))*freq(iharm) + phi(iharm));

		% with magnitude applied
		% temp_pin = temp_pin + A(iharm)*cos(2*pi*t(starti(parti):endi(parti))*freq(iharm) + phi(iharm));
	end
	pin = temp_pin;

	% Generate spectrogram
	window = round(1/decomp_info(parti).F0_actual*60*1000); %round(decomp_info(parti).F0_actual*6); %600;
	ov = window-10; %round(window*0.9833); %590; % or possibly try minus 10 instead?
	[sg,Ftmp,Ttmp] = spectrogram(pin,window,ov,[],fs,'yaxis');
	Ttmp_part(parti) = {Ttmp};
	Ftmp_part(parti) = {Ftmp};
	spect(parti) = {20*log10(abs(sg))};
	
	fi = 1;
	while Ftmp(fi) < max_freq
		fi = fi+1;
	end

end

for parti = 1:num_parts
	this_spect = spect{parti};
	Ttmp = Ttmp_part{parti};
	Ftmp = Ftmp_part{parti};
end


%% Plot

fontsize = 18; 

% Plot magnitude of harmonics 
fig = figure('Renderer', 'painters', 'Position', [100 100 600 365]);
tiledlayout(1, 3, 'TileSpacing','none')

nexttile
stem(freq,20*log10(A),'BaseValue',-inf,'Color','black');
ylabel('Magnitude (dB)')
title('Mag. Spectrum')
xlim([0 max_freq])
ylim([-75 0])
view(90,90)
set(gca, 'XDir','reverse')
set(gca,'FontSize',fontsize)
set(gca,'YTickLabel',[]);
set(gca,'box','off')
hold on
num_harms = 25; 
for iharm = 1:num_harms
	xline(decomp_info(parti).F0_actual*iharm, ':');
end
xticklabels([0 1, 2, 3, 4, 5])
xlabel('Frequency (kHz)')

% Plot 'spectrogram' of phase
ax(2) = nexttile([1, 2]);
pcolor(1000*Ttmp,Ftmp(1:fi),this_spect(1:fi,:));
shading interp
xlabel('Time (ms)')
set(gca,'YTickLabel',[]);
title('Synthesized Spectrogram')
colormap(gray)
ax(2).CLim = [max(max(this_spect(1:fi,:)))-20, max(max(this_spect(1:fi,:)))]; % min was 30
set(gca,'FontSize',fontsize)
set(gca,'box','off')
hold on
xlim([5 1/F0*6*1000])



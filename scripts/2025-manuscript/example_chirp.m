%% example_chirp.m

%% Load in exaple 

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
F0s = pitch_order;
files = files(order);
F0s = F0s(order);

ind = 20; %7;
target_F0 = F0s(ind);
target_file = fullfile(fpath,'waveforms', files{ind});



%% Plot 

% Creates filename for analyzed data
if ismac
	savepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/decomped';
else
	savepath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\decomped';
end
savefile = sprintf('%s_F0_%0.0f_decomp.mat', target, target_F0);

temp = load(fullfile(savepath, savefile));
decomp_info = temp.decomp_info;

% Sets time vector & max frequency
[data,fs] = audioread(target_file);
dur = length(data)/fs; %in sec
t = (1/fs):(1/fs):dur;
t = t(1:length(data)); % not sure why t and data are different lengths, this remedies that

%%%%%% Could change this to be a specific # of harmonics each time 
max_freq = target_F0*25; %3500; % Determines maximum frequency to analyze
if max_freq>fs/2
	max_freq = fs/2-1;
end
%%%%%%

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

		fig = figure('Renderer', 'painters', 'Position', [100 100 1200 800]);
		ax = axes('Position',[0.05 0.56 0.28 0.40]);
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
		set(gca,'FontSize',8)

		ax = axes('Position',[0.38 0.56 0.28 0.40]);
		stem(freq,A);
		xlabel('Frequency (Hz)')
		ylabel('Magnitude')
		title('Approx. Mag. Spectrum')
		xlim([0 max_freq])
		set(gca,'FontSize',8)


		ax = axes('Position',[0.38 0.07 0.28 0.40]);
		stem(freq,uwphi);
		xlabel('Frequency (Hz)')
		ylabel('Phase (radians)')
		title('Approx. Phase Spectrum')
		xlim([0 max_freq])
		set(gca,'FontSize',8)

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

	% Wait on plotting till the very end
	if parti == 1 % spectrogram of original vowel for comparison
		[sg2,Ftmp2,Ttmp2] = spectrogram(data,window,[],[],fs,'yaxis');
		spect2 = 20*log10(abs(sg2));
		ax = axes('Position',[0.71 0.07 0.28 0.40]);
		spec_image2 = pcolor(1000*Ttmp2,Ftmp2(1:fi),spect2(1:fi,:));
		shading interp
		ylabel('Freq (Hz)')
		xlabel('Time (ms)')
		title('Original Spectrogram')
		set(gca,'FontSize',8)
		xlim([2 1/target_F0*10*1000])

	end

end

for parti = 1:num_parts
	this_spect = spect{parti};
	Ttmp = Ttmp_part{parti};
	Ftmp = Ftmp_part{parti};
	ax = axes('Position',[0.71+(0.13825*(parti-1)) 0.56 0.275/num_parts 0.40]);
	spec_image = pcolor(1000*Ttmp,Ftmp(1:fi),this_spect(1:fi,:));
	shading interp
	if parti == 1
		ylabel('Freq (Hz)')
		xlabel('Time (ms)')
		title(['Synthesized Spectrogram, w=' num2str(window)])
		set(gca,'FontSize',8)
		xlim([2 1/target_F0*10*1000])

	else
		set(gca,'YTickLabel',[]);
	end
end

%% Plot 

% Plot magnitude of harmonics
fig = figure('Renderer', 'painters', 'Position', [100 100 600 400]);
ax = axes('Position',[0.12 0.14 0.40 0.80]);
stem(freq,20*log10(A),'BaseValue',-inf,'Color','black');
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title(['Mag. Spectrum, F0=' num2str(round(target_F0)) 'Hz'])
xlim([0 max_freq])
ylim([-75 -15])
view(90,90)
set(gca, 'XDir','reverse')
set(gca,'FontSize',16)
set(gca,'YTickLabel',[]);
set(gca,'box','off')

hold on
num_harms = 25;
for iharm = 1:num_harms
	xline(target_F0*iharm, ':');
end

% Plot 'spectrogram' of phase
ax = axes('Position',[0.54 0.14 0.42 0.80]);
pcolor(1000*Ttmp,Ftmp(1:fi),this_spect(1:fi,:));
shading interp
%ylabel('Frequency (Hz)')
xlabel('Time (ms)')
set(gca,'YTickLabel',[]);
title('Synthesized Spectrogram')

%alternative caxis settings
% colormap(flipud(gray))
% caxis([-40 50])

% Show 10 periods
xlim([10 1/target_F0*10*1000])

colormap(gray)
ax.CLim = [30, max(max(this_spect(1:fi,:)))];
set(gca,'FontSize',16)
set(gca,'box','off')
hold on

%% Export 

savepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/figures/manuscript';
exportgraphics(gcf, fullfile(savepath, 'stimulus-chirp.png'), 'Resolution', 600)


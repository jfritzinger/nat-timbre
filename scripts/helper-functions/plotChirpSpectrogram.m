function plotChirpSpectrogram(vow_str,F0,file)

% Reads in file
[data,fs] = audioread(file);

% Creates filename for analyzed data
if ismac
	savepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/decomped';
else
	savepath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\decomped';
end
savefile = sprintf('%s_F0_%0.0f_decomp.mat', vow_str, F0);

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
MASD = decomp_info(1).MASD{:};

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
		%     freq = (freq(1):freq(1):freq(1)*35)';
		%     phi = zeros(length(phi),1);

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
		altplot(iharm,:) = cos(2*pi*t*freq(iharm) + phi(iharm));

	end
	pin = temp_pin;

	% Generate spectrogram
	% window = round(1/decomp_info(parti).F0_actual*60*1000); %round(decomp_info(parti).F0_actual*6); %600;
	window = round(1/decomp_info(parti).F0_actual*fs); 
	%ov = window-10; %round(window*0.9833); %590; % or possibly try minus 10 instead?
	ov = round(0.9*window);
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
		[sg2,Ftmp2,Ttmp2] = spectrogram(data,window,ov,[],fs,'yaxis');
		spect2 = 20*log10(abs(sg2));
		ax = axes('Position',[0.71 0.07 0.28 0.40]);
		spec_image2 = pcolor(1000*Ttmp2,Ftmp2(1:fi),spect2(1:fi,:));
		shading interp
		ylabel('Freq (Hz)')
		xlabel('Time (ms)')
		title('Original Spectrogram')
		set(gca,'FontSize',8)
		xlim([2 1/F0*10*1000])
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
		xlim([2 1/F0*10*1000])
	else
		set(gca,'YTickLabel',[]);
	end
end
annotation('textbox',[0.1 0 .4 .5],'String',sprintf('%s, %0.0f Hz',vow_str,F0),'FitBoxToText','on')

%% Plot 
fontsize = 8; 

% Plot magnitude of harmonics 
fig = figure('Renderer', 'painters', 'Position', [100 100 600 400]);
ax(1) = axes('Position',[0.12 0.14 0.40 0.80]);
stem(freq,20*log10(A),'BaseValue',-inf,'Color','black');
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('Approx. Mag. Spectrum')
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

% Plot 'spectrogram' of phase 
% ax(2) = axes('Position',[0.54 0.14 0.42 0.80]);
% % spec_image = pcolor(1000*Ttmp,Ftmp(1:fi),this_spect(1:fi,:));
% pcolor(1000*Ttmp,Ftmp(1:fi),this_spect(1:fi,:));
% shading interp
% ylabel('Frequency (Hz)')
% xlabel('Time (ms)')
% set(gca,'YTickLabel',[]);
% title('Synthesized Spectrogram')
% % colormap(flipud(gray)) % alternative caxis settings
% % caxis([-40 50]) % alternative caxis settings
% colormap(gray)
% ax(2).CLim = [max(max(this_spect(1:fi,:)))-20, max(max(this_spect(1:fi,:)))]; % min was 30
% set(gca,'FontSize',fontsize)
% set(gca,'box','off')
% hold on
ax(2) = axes('Position',[0.54 0.14 0.42 0.80]);
f_alt = F0:F0:max_freq;
t_alt = (1/fs):(1/fs):dur;
s = surf(t_alt,f_alt./1000,altplot);
s.EdgeColor = 'none';
view(0,90)
%xlim(xlimits)
title(sprintf('Spectrogram, F0=%0.0f', F0))
xlabel('Time (s)')
ylabel('CFs (kHz)')
%ylim(ylimits)
set(gca, 'FontSize', 16)
shading interp
xlim([0.1 1/F0*2+0.1])

% Calculate MASD and add to plot (JBF)
% MASD = calculateMASD(Ftmp, Ttmp, this_spect, fi, F0, 0);
% f = Ftmp(1:fi);
% ax(3) = axes('Position',[0.83 0.14 0.12 0.80]);
% plot(f(1:end-1), MASD, 'k')
% set(gca, 'XDir','reverse')
% view(90,90)
% hold on
% num_harms = 25;
% for iharm = 1:num_harms
% 	xline(F0*iharm, ':');
% end
% yline(0, 'k')
% xticks([])
% xlim([0 max_freq])
% title('MASD')
% set(gca,'FontSize',fontsize)
% set(gca,'box','off')
%ylim([-1.5 1.5])

% Reposition axes (JBF) 
set(ax(1), 'Position', [0.1  0.14 0.26 0.80])
set(ax(2), 'Position', [0.43 0.14 0.42 0.80])
%set(ax(3), 'Position', [0.83 0.14 0.12 0.80])

end
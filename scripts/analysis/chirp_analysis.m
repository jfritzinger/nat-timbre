function chirp_analysis(vow_str,vow_F0,file)

% Reads in file
[data,fs] = audioread(file);

% Creates filename for analyzed data
if ismac
	savepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/decomped';
else
	savepath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\decomped';
end
savefile = sprintf('%s_F0_%0.0f_decomp.mat', vow_str, vow_F0);

% Sets time vector & max frequency
dur = length(data)/fs; %in sec
t = (1/fs):(1/fs):dur;
t = t(1:length(data)); % not sure why t and data are different lengths, this remedies that

%%%%%% Could change this to be a specific # of harmonics each time 
max_freq = vow_F0*25; %3500; % Determines maximum frequency to analyze
if max_freq>20000
	max_freq = 20000;
end
%%%%%%

opts = optimoptions('fmincon','Display','off');

% Run analysis
if ~isfile(fullfile(savepath, savefile)) % Check if we already have fit this vowel

	%we need to separate the overall vowel into segments so that we can
	%somewhat follow the frequency/phase changes over time
	%first, try separating into halves
	num_parts = 1;
	starti(1) = 1;
	for parti = 1:num_parts
		endi(parti) = ceil((parti*length(data))/num_parts);
		data_parts(parti) = {data(starti(parti):endi(parti))};
		starti(parti+1) = endi(parti) + 1;
	end
	starti(end) = [];

	for parti = 1:num_parts
		data_part = data_parts{parti};
		t_part = t(starti(parti):endi(parti));
		[pxx, f] = pwelch(data_part,[],[],[],fs);
		[pks,locs] = findpeaks(pxx);

		pks_dB = 10*log10(pks);
		dB_thresh = -70;

		i = 0;
		this_pk = -inf;
		while this_pk <= dB_thresh
			i = i + 1;
			this_pk = pks_dB(i);
		end

		F0_actual = f(locs(i));

		harmonics(1) = F0_actual;
		current_F0_guess = F0_actual;

		i = 2;
		this_f = 0;
		while this_f <= max_freq
			next_F_guess = current_F0_guess*i;
			[~,I(i)] = min(abs(next_F_guess - f(locs)));
			if I(i-1)==I(i)
				harmonics(i) = next_F_guess;
				this_f = next_F_guess;
			else
				harmonics(i) = f(locs(I(i)));
				this_f = f(locs(I(i)));
			end
			current_F0_guess = (current_F0_guess*(i-1) + (harmonics(i)-harmonics(i-1)))/i;
			i = i+1;
		end

		% Added in because last harmonic can be above the max_freq in
		% previous while loop
		if harmonics(end)>max_freq
			harmonics(end)=[];
		end

		% Checks that no sporious peaks are being treated as harmonics and
		% gets rid of them if there are any spurious peaks present
		iharm = 1;
		extras = 0;
		while iharm < length(harmonics)
			freqdiff(iharm) = harmonics(iharm+1) - harmonics(iharm);
			if freqdiff(iharm) < harmonics(1)*0.65
				harmonics(iharm+1) = [];
				extras = extras + 1;
				iharm = iharm - 1;
			end
			iharm = iharm + 1;
		end
		fprintf('Removed %i extra components\n',extras)

		% Loop through all harmonics
		for iharm = 1:length(harmonics)
			%generate Kaiser window for FIR filter
			fcuts = [harmonics(iharm)-(harmonics(1)/2) harmonics(iharm)-10 ...
				harmonics(iharm)+10 harmonics(iharm)+(harmonics(1)/2)];
			mags = [0 1 0];
			devs = [0.01 0.05 0.01];

			[n,Wn,beta,ftype] = kaiserord(fcuts,mags,devs,fs);
			n = n + rem(n,2); %odd length means even order
			hh = fir1(n,Wn,ftype,kaiser(n+1,beta),'noscale');

			filtout = filter(hh,1,data_part);
			% figure
			% [H,f] = freqz(hh,1,1024,fs);
			% plot(f,abs(H))
			% grid

			%do "simple method"
			%cross correlation filtout with cosine to find phase.
			%we assume that the frequency aligns with the peak in the power spectral
			%density
			simple_cos = cos(2*pi*t*harmonics(iharm));
			[c,lags] = xcorr(filtout,simple_cos);
			period = 1/harmonics(iharm);


			[m,I] = max(c(length(filtout):length(filtout)+floor((period*fs))));

			delay = lags(length(filtout)+I)/fs; %in s
			simple_phase(iharm) = (delay/period)*pi;

			%         wind_size = 4000;
			%         incr_size = wind_size/2;
			%         start_inds = [1 ceil(length(filtout)/2)];
			%         end_inds = [floor(length(filtout)/2) length(filtout)];
			start_inds = 1;
			end_inds = length(filtout);
			use_whole = 1; %if 1, use whole "filtout", if 0, break into subwindows
			for istart = 1:length(start_inds)
				if use_whole == 1
					this_wind = filtout;
					this_t = t_part;
					end_inds(istart) = length(filtout);
				else
					this_wind = filtout(start_inds(istart):end_inds(istart));
					this_t = t(start_inds(istart):end_inds(istart));
				end
				fun = @(x) find_component(x(1),x(2),x(3),this_wind,this_t);
				lb = [0.00001, 0, -2*pi];
				ub = [max(data_part),max_freq+harmonics(1),2*pi];
				fitreps = 300;
				for ifit = 1:fitreps
					rng('shuffle')
					x0 = ((ub-lb).*rand(1,3))+lb;
					x0(2) = harmonics(iharm);
					[temp_out(ifit,:),temp_fval(ifit),~,~,~,~,~] = fmincon(fun,x0,[],[],[],[],lb,ub,[],opts);
				end
				[~,min_i] = min(temp_fval);
				out = temp_out(min_i,:);
				fval = temp_fval(min_i);
				A(iharm,istart) = out(1);
				freq(iharm,istart) = out(2);
				if out(3) < -pi
					out(3) = out(3) + 2*pi;
				elseif out(3) > pi
					out(3) = out(3) - 2*pi;
				end
				phi(iharm,istart) = out(3);
				all_fvals(iharm,istart) = temp_fval(min_i);
			end
			fprintf('Found fit for component %i of %i\n',iharm,length(harmonics))
		end

		decomp_info(parti).A = {A};
		decomp_info(parti).freq = {freq};
		decomp_info(parti).phi = {phi};
		decomp_info(parti).fitreps = fitreps;
		decomp_info(parti).fval = {all_fvals};
		decomp_info(parti).max_freq = max_freq;
		decomp_info(parti).start_inds = {starti};
		decomp_info(parti).end_inds = {endi};
		decomp_info(parti).num_parts = num_parts;
		decomp_info(parti).F0_actual = F0_actual;

	end

	save(fullfile(savepath, savefile),'decomp_info')
else
	temp = load(fullfile(savepath, savefile));
	decomp_info = temp.decomp_info;
end

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

		%     freq = (freq(1):freq(1):freq(1)*35)';
		%     phi = zeros(length(phi),1);

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
		xlim([2 1/vow_F0*10*1000])

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
		xlim([2 1/vow_F0*10*1000])

	else
		set(gca,'YTickLabel',[]);
	end
end


% spect_parts{parti} = spect;
% Ttmp_parts{parti} = Ttmp2;


% figure
% merged_spect = [spect_parts{1} spect_parts{2}];
% full_Ttmp = [Ttmp_parts{1} Ttmp_parts{2}+Ttmp_parts{1}(end)];
% spec_image = pcolor(1000*full_Ttmp,Ftmp(1:fi),merged_spect(1:fi,:));
% shading interp
% ylabel('Freq (Hz)')
% xlabel('Time (ms)')
% title('Spectrogram')
annotation('textbox',[0.1 0 .4 .5],'String',sprintf('%s, %0.0f Hz',vow_str,vow_F0),'FitBoxToText','on')

% Plot magnitude of harmonics 
fig = figure('Renderer', 'painters', 'Position', [100 100 600 400]);
ax = axes('Position',[0.12 0.14 0.40 0.80]);
stem(freq,20*log10(A),'BaseValue',-inf,'Color','black');
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('Approx. Mag. Spectrum')
xlim([0 max_freq])
ylim([-75 0])
view(90,90)
set(gca, 'XDir','reverse')
set(gca,'FontSize',8)
set(gca,'YTickLabel',[]);
set(gca,'box','off')

hold on
num_harms = 25; 
for iharm = 1:num_harms
	xline(decomp_info(parti).F0_actual*iharm, ':');
end

% Plot 'spectrogram' of phase 
ax = axes('Position',[0.54 0.14 0.42 0.80]);
% spec_image = pcolor(1000*Ttmp,Ftmp(1:fi),this_spect(1:fi,:));
pcolor(1000*Ttmp,Ftmp(1:fi),this_spect(1:fi,:));
shading interp
ylabel('Frequency (Hz)')
xlabel('Time (ms)')
set(gca,'YTickLabel',[]);
title('Synthesized Spectrogram')

%alternative caxis settings
% colormap(flipud(gray))
% caxis([-40 50])

% Show 10 periods 
xlim([2 1/vow_F0*10*1000])

colormap(gray)
ax.CLim = [max(max(this_spect(1:fi,:)))-20, max(max(this_spect(1:fi,:)))]; % min was 30
set(gca,'FontSize',8)
set(gca,'box','off')
hold on

end

%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sse = find_component(A,f,phi,filtout,t)
thisfit = (A*cos(2*pi*t*f + phi));

sse = sum((thisfit - filtout').^2);%sum of squared errors
end
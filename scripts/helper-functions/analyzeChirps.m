function analyzeChirps(vow_str,F0,file)

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

opts = optimoptions('fmincon','Display','off');

% Run analysis
if F0~=0%~isfile(fullfile(savepath, savefile)) % Check if we already have fit this vowel

	%we need to separate the overall vowel into segments so that we can
	%somewhat follow the frequency/phase changes over time
	%first, try separating into halves
	num_parts = 1;
	endi = NaN(num_parts,1);
	starti = NaN(num_parts+1, 1);
	starti(1) = 1;
	for parti = 1:num_parts
		endi(parti) = ceil((parti*length(data))/num_parts);
		data_parts(parti) = {data(starti(parti):endi(parti))};
		starti(parti+1) = endi(parti) + 1;
	end
	starti(end) = [];

	% Initialize struct
	default_struct = struct('A', {}, 'freq', {}, 'phi', {}, 'fitreps', {}, 'fval', {}, ...
		'max_freq', {}, 'start_inds', {}, 'end_inds', {}, ...
		'num_parts', {}, 'F0_actual', {}, 'MASD', {});
	decomp_info = repmat(default_struct, 1, num_parts);
	for parti = 1:num_parts
		data_part = data_parts{parti};
		t_part = t(starti(parti):endi(parti));
		[pxx, f] = pwelch(data_part,[],[],[],fs);
		[pks,locs] = findpeaks(pxx);

		pks_dB = 10*log10(pks);
		dB_thresh = -70;
		i = find(pks_dB > dB_thresh, 1);
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
		% JBF updated to remove while loop 
		n_harmonics = length(harmonics);
		to_remove = zeros(1, n_harmonics);
		freqdiff = diff(harmonics);
		to_remove(2:end) = freqdiff < harmonics(1) * 0.65;
		harmonics = harmonics(~to_remove);
		extras = sum(to_remove);
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
				parfor ifit = 1:fitreps
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

		% construct "not-spectrogram" spectrogram using phase first, append
		% waveforms of cosines with appropriate phase into the same array
		phase_spect = zeros(length(decomp_info.freq{1}),length(t));
		for icomp = 1:length(decomp_info.freq{1})
			A2 = decomp_info.A{1}(icomp);
			f2 = decomp_info.freq{1}(icomp);
			phi2 = decomp_info.phi{1}(icomp);
			phase_spect(icomp,:) = (A2*cos(2*pi*t*f2 + phi2));
		end

		% Manually unwrap phase to be monotonically decreasing
		% JBF, edited to get rid of while loops (I just don't like reading
		% them) 
		if parti == 1
			A = decomp_info(parti).A{:};
			freq = decomp_info(parti).freq{:};
			phi = decomp_info(parti).phi{:};
			uwphi2 = zeros(size(phi));
			uwphi2(1) = phi(1);
			for i = 2:length(phi)
				uwphi2(i) = phi(i) - 2*pi * ceil((phi(i) - uwphi2(i-1)) / (2*pi));
			end
		end

		% Reconstruct signal using these phases
		temp_pin = zeros(1,length(t(starti(parti):endi(parti))));
		for iharm = 1:length(freq)
			temp_pin = temp_pin + cos(2*pi*t(starti(parti):endi(parti))*...
				freq(iharm) + phi(iharm)); % without magnitude		
			% temp_pin = temp_pin + A(iharm)*cos(2*pi*t(starti(parti):...
			% 	endi(parti))*freq(iharm) + phi(iharm)); % with magnitude applied
		end
		pin = temp_pin;

		% Generate spectrogram
		window = round(1/decomp_info(parti).F0_actual*60*1000); %round(decomp_info(parti).F0_actual*6); %600;
		ov = window-10; %round(window*0.9833); %590; % or possibly try minus 10 instead?
		[sg,Ftmp,Ttmp] = spectrogram(pin,window,ov,[],fs,'yaxis');
		spect(parti) = {20*log10(abs(sg))};
		fi = find(Ftmp >= max_freq, 1);

		% Calculate MASD and add to plot
		MASD = calculateMASD(Ftmp, Ttmp, spect{parti}, fi, F0, 0);
		decomp_info(parti).MASD = {MASD};
	end
	save(fullfile(savepath, savefile),'decomp_info')
else
	disp('Analysis already saved!')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sse = find_component(A,f,phi,filtout,t)
thisfit = (A*cos(2*pi*t*f + phi));

sse = sum((thisfit - filtout').^2);%sum of squared errors
end
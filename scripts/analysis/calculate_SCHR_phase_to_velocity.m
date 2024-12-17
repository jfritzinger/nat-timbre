%% calculate_SCHR_phase_to_velocity
clear

%% Create SCHR

params.type = 'SCHR';
params.version = 1;
params.dur = 0.4;
params.reptim = 1;
params.ramp_dur = 0.025;
params.binmode = 2;
params.stimdB = 60;
params.F0 = 50;
%params.F0 = [25, 50, 100, 200, 400, 600];
%params.C = [-1  0 1];
params.C = -1;
params.nrep = 30;
params.Fs = 100000;
params.mnrep = 1;
params = generate_SCHR(params);

%% Plot SCHR

figure('Position',[54,453,1513,420])
tiledlayout(1, 3)
stim = params.stim;

% Time
nexttile
t = linspace(0, params.dur, params.dur*params.Fs);
plot(t, stim)
xlim([0.02 0.07])
xlabel('Time (s)')

% FFT
nexttile
y2 = fft(stim);
m = abs(y2);
mdB = 20*log10(m);
freq = (0:length(y2)-1)*params.Fs/length(y2);
mdB(mdB<0) = 0;
freq(freq>params.Fs/2) = [];
mdB = mdB(1:length(freq))';
plot(freq, mdB)
xlim([0 1000])
xlabel('Frequency (Hz)')
ylabel('Amplitude (dB SPL)')

% Spectrogram
nexttile
spec = spectrogram(stim, 200, [], freq, params.Fs);
spec = 10*log10(abs(spec));
imagesc(t, freq, spec)
xlabel('Time (s)')
set(gca, 'YDir', 'normal');
ylim([0 16000])
xlim([0 0.1])
ylabel('Frequency (Hz)')

% Plot Phase
% y2 = fft(stim);
% m = abs(y2);
% mdB = 20*log10(m);
% freq = (0:length(y2)-1)*params.Fs/length(y2);
% mdB(mdB<0) = 0;
% freq(freq>params.Fs/2) = [];
% mdB = mdB(1:length(freq))';
%
% theta = angle(y2);
% theta = theta(1:length(theta)/2+1);
%
% figure;
% plot(freq, unwrap(theta));
% xlabel('Frequency (Hz)');
% ylabel('Phase (radians)');
% title('Phase Spectrum');
% grid on;

%% Convert into velocity numerically
% Answer: 0.80 kHz/ms

f_max = 16000;
frequency_span = (f_max - params.F0)/1000;
dur = -1*1/params.F0*params.C*1000;
velocity = frequency_span / dur;

%% Convert into velocity from spectrogram

fs = 100000;
[s, f, t] = spectrogram(stim, hamming(256), 128, 512, fs);
s = 10*log10(abs(s));
s(s<0) = 0;
[ifq, t_if] = instfreq(s, f, t);
t_if = t_if*1000;
ifq = ifq/1000;

dt = diff(t_if);
difq = diff(ifq); %, 1, 2); % Difference along the time axis
velocity = difq ./ dt; % Velocity in Hz/s


% Plot instantaneous frequency
figure('Position',[560,495,950,353])
tiledlayout(1, 2)
nexttile
plot(t_if, ifq)
title('Instantaneous Frequency')

% Plot velocity vs time
nexttile
plot(t_if(2:end), velocity)
hold on
yline(median(velocity, 'omitnan'), '--r', 'LineWidth',2)
xlabel('Time (ms)');
ylabel('Velocity (kHz/ms)');
title('Frequency Sweep Velocity')

%% Use threshold to find velocity
% 
% figure
% fs = 100000;
% [s, f, t] = spectrogram(stim, hamming(256), 128, 512, fs);
% s = 10*log10(abs(s));
% %s(s<0) = 0;
% 
% imagesc(t, freq, s)
% xlabel('Time (s)')
% set(gca, 'YDir', 'normal');
% ylim([0 16000])
% xlim([0 0.1])
% ylabel('Frequency (Hz)')
% %vel = diff(spec);
% 
% max_spec = find(max(s));

%% Find velocity frequency by frequency

f_array = 0:2000:16000;
for isteps = 2:8

	[s, f, t] = spectrogram(stim, 200, [], 512, fs);
	s = 10*log10(abs(s));

	freq_lo = f_array(isteps-1);
	freq_hi = f_array(isteps);
	f_ind = find(f<freq_hi & f>freq_lo);

	s_part = s(f_ind,:);
	period = 1/params.F0;
	period_array = 0:period:params.dur;

	% figure
	% imagesc(t, f(f_ind), s_part)
	% hold on
	% for ii = 1:length(period_array)
	% 	xline(period_array(ii), 'LineWidth',2)
	% end
	% set(gca, 'YDir', 'normal');


	t_ind = 1:length(t); %find(t<period_array(3) & t>period_array(2));
	s_new = s_part(:,t_ind);
	t_new = t(t_ind);
	f_new = f(f_ind);
	%s_new(s_new < 0) = 0;

	% figure
	% imagesc(t_new, f_new, s_new)
	% hold on
	% set(gca, 'YDir', 'normal');

	s_new(s_new<0) = 0;
	[ifq, t_if] = instfreq(s_new, f_new, t_new);
	t_if = t_if*1000;
	ifq = ifq/1000;
	% figure
	% plot(t_if, ifq)

	dt = diff(t_if);
	difq = diff(ifq); %, 1, 2); % Difference along the time axis
	velocity = difq ./ dt; % Velocity in Hz/s
	vel(isteps) = max(velocity);

	% figure
	% plot(velocity)

end

figure
plot(f_array(2:end), vel)

%% Other options

% [~, maxIdx] = max(abs(s_new), [], 1);
% ifq = f_new(maxIdx);
% t_if = t_new'*1000;
% ifq = ifq/1000;
% figure
% plot(t_new, ifq)
%
% z = hilbert(s_new);
% instfrq = fs/(2*pi)*diff(unwrap(angle(z)));
%
% figure
% plot(t_new(2:end),instfrq)
% ylim([0 fs/2])
%
% fridge = tfridge(s_new,f_new);
%
%
% % Average each period together
% % Calculate the period
% T = 1/params.F0; % seconds
% samples_per_period = round(T * (length(t_new) / (t_new(end) - t_new(1))));
% num_periods = floor(length(ifq) / samples_per_period);
%
% folded_psth = reshape(ifq(1:num_periods*samples_per_period), samples_per_period, num_periods);
% period_psth = mean(folded_psth, 2);
% period_time = linspace(0, T, samples_per_period);
%
% % Plot the period PSTH
% figure;
% plot(period_time, period_psth);
% xlabel('Time within period (s)');
% ylabel('Spike rate');
% title('Period PSTH (f0 = 100 Hz)');
%
% dt = diff(t_if);
% difq = diff(ifq); %, 1, 2); % Difference along the time axis
% velocity = difq ./ dt; % Velocity in Hz/s



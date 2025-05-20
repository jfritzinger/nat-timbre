%% Calculating ERR and modulation depth of instrument sounds 

%% Calculate ERR 

% Simulate a modulated signal
fs = 1000; % Sampling rate (Hz)
t = 0:1/fs:1; % Time vector
carrier_freq = 100; % Hz
modulation_freq = 10; % Hz (ERR)
signal = sin(2*pi*carrier_freq*t) .* (1 + 0.5*sin(2*pi*modulation_freq*t));
envelope = abs(hilbert(signal)); % MATLAB example


% Compute envelope spectrum
[es, f] = envspectrum(signal, fs, 'Band', [1 100]); 
%[es, f] = envspectrum(signal, fs); % MATLAB's envspectrum function [4]


% Find dominant ERR
[~, idx] = max(es(2:end)); % Skip DC (0 Hz)
ERR = f(idx + 1); % e.g., 10 Hz



%% Calculate modulation depth 


% Example signal (replace with your own)
t = 0:1e-4:0.1;
x = (1 + 0.5*cos(2*pi*50*t)) .* cos(2*pi*1000*t); % 50 Hz modulation, 1000 Hz carrier

% Extract envelope using Hilbert transform
env = abs(hilbert(x)); % [1][2]
% OR
[env, ~] = envelope(x, 'analytic'); % [2]

% Modulation depth 
Emax = max(env);
Emin = min(env);
modDepth = (Emax - Emin) / (Emax + Emin);

% Plot envelope 
plot(t, x); hold on;
plot(t, env, 'r', 'LineWidth', 2);
legend('Signal', 'Envelope');
hold off;

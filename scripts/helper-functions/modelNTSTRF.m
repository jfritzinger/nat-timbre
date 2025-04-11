function [R2, avModel, stdModel, ratio, max_all, R] = modelNTSTRF(param, data_STRF, data)
% Function that uses STRF to predict responses to WBTIN
% J. Fritzinger, created 7/10/23

% Get STRF data 
STRF_mat = data_STRF.H2ex_strf-data_STRF.H2in_strf;
f = data_STRF.f;
f(f>16000) = [];
STRF_mat(length(f)+1:end, :) = [];

% Recreate stimulus
param.physio = 1;
if ~isfield(param, 'target')
	param.target = extractBefore(param.filename{1}, '.');
end
[param] = generate_NT(param);

% Calculate the spectrogram of a stimulus
windowLength = 200; %50; % Length of the analysis window (in samples)
noverlap = 199; %49; % Overlap size between consecutive windows (in samples)
nfft = 960; % Length of the FFT

% Convolve STRF and the stimulus
num_stim = length(param.list);
all_avg_rate = zeros(1, num_stim);
for istim = 1:num_stim

	% Plot the spectrogram of a stimulus
	[S,freq, ~] = spectrogram(param.stim(istim, :), windowLength, noverlap, nfft, param.Fs);
	freq(freq>16000) = [];
	S(length(freq)+1:end, :) = [];
	S = abs(S);

	% Plot 

	% Convolve
	for ifreq = 1:length(freq)
		signal = S(ifreq,:);
		filter = STRF_mat(ifreq,:);
		conv1 = conv(signal, filter, 'same');
		convolution_result(ifreq,:) = conv1;
	end

	% Get PSTH & average rate
	psth = sum(convolution_result, 1);
	avg_rate = mean(psth);
	all_avg_rate(istim) = avg_rate;
end

[Timbres,~,Timbrei] = unique([param.list.iwav_file].');
num_Timbres = length(Timbres);
rate_size = [num_Timbres,1]; % [num_F0s,num_Fps];
[avModel,stdModel,~,~] = accumstats({Timbrei},all_avg_rate, rate_size);

% Model
max_rate = max(data.rate);
max_rate_m = max(avModel);
ratio =  max_rate/max_rate_m;
max_all = max([max_rate, max((avModel.*ratio))]);

% Measure prediction accuracy
R_int = corrcoef(avModel,data.rate);
R2 = R_int(1, 2).^2;
R = R_int(1,2);

end

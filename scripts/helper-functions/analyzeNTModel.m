function data = analyzeNTModel(model, param, MTF_type, CF)
% plotNT plots the model response to the natural timbre stimuli. 
%
% Inputs:
%
% Outputs:
%
%
% Other m-files required: accumstats
%
% Depends on 'Tuning.xlsx' excel file
%
% Author: J. Fritzinger
% Created: 2023-09-27; Last revision: 2024-09-16
%
% -------------------------------------------------------------------------

% Get BE or BS model outputs
if strcmp(MTF_type, 'BE')
	spike_rates = model.average_ic_sout_BE;
	model_response = squeeze(model.ic_BE);

elseif strcmp(MTF_type, 'BS')
	spike_rates = model.average_ic_sout_BS;
	model_response = squeeze(model.ic_BS);

end

% Downsample
dur = 300;
fs = 100000;
fs_new = 4000;
D = fs / fs_new;
model_response = resample(model_response', 1, D)';
nsamples_new = size(model_response, 2);
t = (0:nsamples_new-1) / fs_new*1000;
fs = fs_new;
ind_stim = t<dur;
model_response = model_response(:,ind_stim);
t = t(ind_stim);


% Load in tuning sheet
base = getPathsNT();
tuning = readtable(fullfile(base, 'Tuning.xlsx'));


% Rate analysis
[Timbres,~,Timbrei] = unique([param.mlist.iwav_file].');
num_Timbres = length(Timbres);
rate_size = [num_Timbres,1]; % [num_F0s,num_Fps];
[rate,rate_std] = accumstats({Timbrei},spike_rates, rate_size); % uses onset?


% Sort by frequency of pitch
index = NaN(num_Timbres, 1);
note_names = extractBetween(param.filename, 'ff.','.wav');
for ii = 1:num_Timbres % Find index of each note in tuning spreadsheet
	if length(note_names{ii}) > 3
		note_names{ii} = extractAfter(note_names{ii}, '.');
	end
	index(ii) = find(strcmp(note_names(ii), tuning.Note));
end
pitch_order = tuning.Frequency(index); % Get freqs of each note
[~, order] = sort(pitch_order); % Sort freqs
pitch = categorical(round(pitch_order(order))); % Order pitches
rate = rate(order);
if isfield(param, 'target')
	F0s = calcManualF0s(param.target);
else
	if num_Timbres==40
		F0s = calcManualF0s('Bassoon');
	elseif num_Timbres==35
		F0s = calcManualF0s('Oboe');
	end
end

% Calculate average rate for each repetition and then sort
reps = length(Timbrei)/num_Timbres;
raw_rates = zeros(reps,num_Timbres);
for j = 1:num_Timbres
	raw_rates(:,j) = spike_rates(j==Timbrei);
end
raw_rates = raw_rates(:,order);

% Get spike time information in similar format as data 
reps = zeros(num_Timbres,1);
for j1 = 1:num_Timbres
	j = Timbrei(:,1) == order(j1);
	reps(j) = (1:sum(j)).';
end

PSTH = zeros(num_Timbres, size(model_response, 2));
for j1 = 1:num_Timbres
	k = Timbrei == order(j1);
	PSTH_allreps{j1} = model_response(k,:);
	PSTH(j1,:) = mean(model_response(k,:), 1);
end



% Saved data
data.rate = rate; 
data.rate_std = rate_std(order);
data.F0s = pitch_order;
data.note_names = note_names;
data.order = order;
data.pitch = pitch;
data.pitch_num = pitch_order(order);
data.raw_rates = raw_rates;
data.F0s_actual = F0s;
data.PSTH_all_reps = PSTH_allreps;
data.PSTH = PSTH;
data.t = t;

%% Temporal Analysis

for i_stim = 1:num_Timbres

	% Parameters 
	dur = 300; % duration in ms

	% Cut off onset from PSTH
	onset = 25; % 25 ms onset 
	ind_onset = t<onset;
	spike_times = PSTH(i_stim,~ind_onset);
	t_times = t(~ind_onset)-onset;

	% Truncate spikes that don't occur in a full cycle of the stimulus
	period = 1000 / F0s(i_stim);	% Period in ms
	num_periods = floor((dur-onset)/period);	% Number of full periods in the stimulus
	samples_per_period = floor(fs*period/1000);
	subset_spike_times = spike_times(1:samples_per_period*num_periods);
	t_subset = t_times(1:samples_per_period*num_periods);

	% Plots
	% figure
	% plot(t_subset, subset_spike_times)
	% hold on
	% periods = 0:period:275;
	% xline(periods)

	% Period PSTH
	period_hist = reshape(subset_spike_times, samples_per_period,[]);
	counts = mean(period_hist, 2);
	edges = linspace(0, period, fs*period/1000);

	% Calculate vector strength
	f = F0s(i_stim);
	VS = abs(1/sum(subset_spike_times) * sum(subset_spike_times .* exp(1i * 2*pi * f .* t_subset)));

	% Calculate VS to CF
	VS_CF = abs(1/sum(subset_spike_times) * sum(subset_spike_times .* exp(1i * 2*pi * CF .* t_subset)));

	% Calculate reliability metric (with onset)
	PSTH_all = PSTH_allreps{i_stim};
	PSTH_even = mean(PSTH_all(2:2:20,:));
	PSTH_odd = mean(PSTH_all(1:2:20,:));
	r = corrcoef(PSTH_even, PSTH_odd);
	r_splithalf = r(1,2);

	% Save outputs
	data.p_hist{i_stim,:} = counts;
	data.t_hist{i_stim,:} = edges;
	data.VS(i_stim) = VS;
	data.VS_CF(i_stim) = VS_CF;
	data.r_splithalf(i_stim) = r_splithalf;
end

end


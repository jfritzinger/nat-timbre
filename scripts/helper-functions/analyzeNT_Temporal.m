function temporal = analyzeNT_Temporal(data_ST, CF)

num_stim = length(data_ST.pitch_num);
for i_stim = 1:num_stim

	% Calculate rasters
	x = data_ST.spike_times{i_stim};
	y = data_ST.spike_reps{i_stim};
	valid = x<0.3e6;
	x = x(valid);
	y = y(valid);

	% PSTH and smoothed PSTH
	edges = linspace(0, 300000,501);
	[PSTH, t] = histcounts(x, edges);
	PSTH_smooth = smooth(PSTH);

	% Calculate period histogram
	% Cut off onset 
	freq = data_ST.pitch_num(i_stim);
	period = 1000 / freq; % Period in ms
	num_bins = 30; % Number of bins for histogram
	wrapped_times = mod(x/1000, period); % Wrap spike times to one period
	edges = linspace(0, period, num_bins+1); % Bin edges
	counts = histcounts(wrapped_times, edges); % Create histogram

	% Normalize counts to get firing rate (spikes/sec)
	% bin_width = period / num_bins; % Bin width in ms
	% firing_rate = counts / (30 * bin_width / 1000); % Normalize by trials and bin width

	% Example spike times (300 ms window, 200 Hz stimulus)
	spike_times = x/1000;
	freq = data_ST.pitch_num(i_stim);

	% Exclude 50ms onset
	spike_times(spike_times<50) = [];

	% Calculate vector strength
	period = 1000 / freq;
	phases = 2 * pi * mod(spike_times, period) / period;
	VS = abs(mean(exp(1i * phases)));
	%vectors = exp(-1i*2*pi*200*spike_times/1000); % same equation
	%sync(j) = abs(mean(vectors)); % same equation

	period = 1000 / CF;
	phases = 2 * pi * mod(spike_times, period) / period;
	VS_CF = abs(mean(exp(1i * phases)));

	% Calculate vector strength for each repetition of the stimulus
	period = 1000 / freq;
	nreps = 20;
	for ind = 1:nreps
		rep_ind = y==ind;
		raw_spikes = x(rep_ind)/1000;
		phases = 2 * pi * mod(raw_spikes, period) / period;
		VS2(ind) = abs(mean(exp(1i * phases)));
		% figure
		% nexttile
		% plot(phases)
		% title(num2str(VS2(ind)))
		% nexttile
		% histogram(raw_spikes, 301)
	end

	% Save outputs
	temporal.x{i_stim} = x;
	temporal.y{i_stim} = y;
	temporal.t = t/1000;
	temporal.PSTH(i_stim,:) = PSTH;
	temporal.PSTH_smooth(i_stim,:) = PSTH_smooth;
	temporal.p_hist(i_stim,:) = counts;
	temporal.t_hist(i_stim,:) = edges;
	temporal.VS(i_stim) = VS;
	temporal.VS_CF(i_stim) = VS_CF;
	temporal.VS_rep(i_stim, :) = VS2;
end

end
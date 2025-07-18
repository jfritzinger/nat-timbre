function temporal = analyzeNT_Temporal(data_ST, CF)
% figure
% tiledlayout(40, 1, "TileSpacing","none")
num_stim = length(data_ST.pitch_num);
for i_stim = 1:num_stim

	% Calculate rasters
	dur = 300; % duration in ms
	nrep = 20; % number of stimulus repetitions
	x = data_ST.spike_times{i_stim};
	y = data_ST.spike_reps{i_stim};
	valid = x<0.3e6;
	x = x(valid)/1000;
	y = y(valid);

	% Calculate spike times per rep
	spike_times = cell(20, 1);
	for ii = 1:20
		spike_times{ii} = x(y==ii)';
	end

	% PSTH (0.25 ms bins)
	num_bins = 1200;
	edges_psth = linspace(0, dur,num_bins+1);
	[PSTH, t] = histcounts(x, edges_psth);

	% Get spike times without onset
	onset = 25; % 25 ms onset
	ind_onset = x<onset;
	spike_times_on = x(~ind_onset);
	reps = y(~ind_onset);

	% Truncate spikes that don't occur in a full cycle of the stimulus
	period = 1000 / (data_ST.F0s_actual(i_stim));	% Period in ms
	num_periods = floor((dur-onset)/period);	% Number of full periods in the stimulus
	ind_full_period = spike_times_on<(num_periods*period);
	subset_spike_times = spike_times_on(ind_full_period);
	subset_reps = reps(ind_full_period);

	% Calculate period histogram
	% Try 0.25 ms bins
	wrapped_times = mod(subset_spike_times, period);
	num_bins = 30; %round(period/0.25); % 30; % Number of bins for histogram
	edges = linspace(0, period, num_bins+1); % Bin edges
	counts = histcounts(wrapped_times, edges); % Create histogram

	% % Normalize counts to get firing rate (spikes/sec)
	% bin_width = period / num_bins; % Bin width in ms
	% firing_rate = counts / (nrep * bin_width); % Normalize by trials and bin width

	% % Plots
	% figure
	% scatter(x, y)
	% hold on
	% periods = 25:period:300;
	% xline(periods)

	% Calculate vector strength
	phases = 2 * pi * (mod(subset_spike_times, period) / period);
	VS = abs(mean(exp(1i * phases)));
	if ~isempty(phases)
		p_value = circ_rtest(phases); % Rayleigh statistic (P < 0.01)
	else
		p_value = NaN;
	end

	% Calculate vector strength for each repetition of the stimulus
	VS_perrep = NaN(nrep,1);
	p_value_perrep = NaN(nrep,1);
	for ind = 1:nrep
		raw_spikes = subset_spike_times(subset_reps==ind);
		phases = 2 * pi * (mod(raw_spikes, period) / period);
		VS_perrep(ind) = abs(mean(exp(1i * phases)));
		if ~isempty(phases)
			p_value_perrep(ind) = circ_rtest(phases); % Rayleigh statistic (P < 0.01)
		end
	end

	% Calculate vector strength for 1-6 harmonics
	for iharm = 1:30
		harm = (data_ST.F0s_actual(i_stim)+(data_ST.F0s_actual(i_stim)*(iharm-1)));
		period = 1000 / harm; % Get period of each harmonic
		phases = 2 * pi * (mod(subset_spike_times, period) / period);
		VS_harms(iharm) = abs(mean(exp(1i * phases)));
		if ~isempty(phases)
			p_value_harms(iharm) = circ_rtest(phases); % Rayleigh statistic (P < 0.01)
		else
			p_value_harms(iharm) = NaN;
		end
		harms(iharm) = harm;
	end

	% Find max phase locking frequency
	% VS_freq_max = 0;
	% all_freqs = linspace(1, 11760, 11760*2);
	% for i = 1:11760*2
	% 	freq = all_freqs(i);
	% 	period = 1000 / freq; % Get period of each harmonic
	% 	phases = 2 * pi * (mod(subset_spike_times, period) / period);
	% 	VS = abs(mean(exp(1i * phases)));
	% 	VS_freq_max(i) = VS;
	% end
	% nexttile
	% hold on
	% % for iharm = 1:10
	% % 	harm = (data_ST.F0s_actual(i_stim)+(data_ST.F0s_actual(i_stim)*(iharm-1)));
	% % 	xline(harm, 'r')
	% % end
	% plot(all_freqs, VS_freq_max)
	% xlim([1 800])

	% Calculate vector strength for 1-6 harmonics
	% for iharm = 1:30
	% 	harm = round((data_ST.F0s_actual(i_stim)+(data_ST.F0s_actual(i_stim)*(iharm-1))));
	% 
	% 	% Get range 
	% 	harm_range = [harm-10 harm+10];
	% 	[~, lo_ind] = min(abs(all_freqs-harm_range(1)));
	% 	[~, hi_ind] = min(abs(all_freqs-harm_range(2)));
	% 	VS_range = VS_freq_max(lo_ind:hi_ind);
	% 	max_VS = max(VS_range);
	% 	VS_harms(iharm) = max_VS;
	% end


	% Calculate ISI
	% isi_all = [];
	% for ind = 1:nrep
	% 	isi = diff(spike_times(reps==ind));
	% 	isi_all = [isi_all; isi]; % will have 20 less indices than spike_times
	% end
	ISI = arrayfun(@(ii) diff(spike_times_on(reps==ii)), 1:nrep, 'UniformOutput', false);
	ISI_all = vertcat(ISI{:});

	% Calculate ISI histogram
	nbins = 80; % for 0.25 ms resolution
	edges_isi = linspace(0, 20, nbins+1);
	ISI_counts_all = histcounts(ISI_all, edges_isi);

	% Calculate VS to CF
	period_CF = 1000 / CF;
	phases = 2 * pi * (mod(spike_times_on, period_CF) / period_CF);
	VS_CF = abs(mean(exp(1i * phases)));

	% Calculate reliability metric (with onset)
	num_bins = 1200;
	edges2 = linspace(0, dur, num_bins+1);
	rep_ind_even = mod(y,2) == 0;
	raw_spikes_even = x(rep_ind_even);
	[PSTH_even, ~] = histcounts(raw_spikes_even, edges2);

	rep_ind_odd = mod(y,2) == 1;
	raw_spikes_odd = x(rep_ind_odd);
	[PSTH_odd, ~] = histcounts(raw_spikes_odd, edges2);
	r = corrcoef(PSTH_even, PSTH_odd);
	r_splithalf = r(1,2);


	% Save outputs
	temporal.x{i_stim} = x;
	temporal.y{i_stim} = y;
	temporal.t = t/1000;
	temporal.PSTH(i_stim,:) = PSTH;
	temporal.p_hist{i_stim,:} = counts;
	temporal.t_hist{i_stim,:} = edges;
	temporal.VS(i_stim) = VS;
	temporal.VS_p(i_stim) = p_value;
	temporal.VS_CF(i_stim) = VS_CF;
	temporal.VS_perrep(i_stim, :) = VS_perrep;
	temporal.VS_p_perrep(i_stim, :) = p_value_perrep;
	temporal.r_splithalf(i_stim) = r_splithalf;
	temporal.ISI_all{i_stim} = ISI_all;
	temporal.ISI_counts_all(i_stim,:) = ISI_counts_all;
	temporal.ISI_edges = edges_isi;
	temporal.VS_harms(i_stim,:) = VS_harms;
	temporal.VS_p_harmns(i_stim,:) = p_value_harms;
	temporal.harms(i_stim,:) = harms;
	temporal.spike_times{i_stim} = spike_times;

end

if num_stim == 0
	temporal = 0;
end

end
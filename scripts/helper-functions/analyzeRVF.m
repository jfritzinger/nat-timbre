function data = analyzeRVF(params)

% Grab one (singluar) set of values for various parameters out of (plural) cell arrays for this particular DSID
param = params;
stimx = param.stim;
cluster = param.cluster;


fs_RVF = 48000; %this is saved from metafile in post_process
jitter_range = param.jitter_range;
if isfield(param,'F0')
	F0s = param.F0;
elseif isfield(param,'F0s')
	F0s = param.F0s;
end
num_F0s = length(F0s);
if isfield(param,'C')
	Cs = param.C;
elseif isfield(param,'Cs')
	Cs = param.Cs;
end
high_freq_limit = param.high_freq_limit;
SPL = param.spl;
nrep = param.nrep;
seg_reps = param.seg_reps;
per_reps = param.per_reps;

%generate a random order of chirp velocities for each burst rep.
%also, generate a length of jitter times
rng(param.seed.Seed) %use stimulus seed
num_vels = length(F0s) * length(Cs);
for i = 1:param.seg_reps
	vel_list(((i-1)*num_vels)+1:i*num_vels) = randperm(num_vels);
	jitter_list(((i-1)*num_vels)+1:i*num_vels) = round(((ceil(rand(1,num_vels)*(jitter_range(2)-jitter_range(1)))+jitter_range(1))/1000)*fs_RVF);
end

%generate the length of all schr waveforms
for F0i = 1:length(F0s)
	for Ci = 1:length(Cs)
		chirp_dur = param.per_reps/F0s(F0i); %duration of the current schroeder period AKA per_reps/frequency
		chirp_npts(F0i+length(F0s)*(Ci-1)) = floor(chirp_dur .* fs_RVF);
		velocities(F0i+length(F0s)*(Ci-1)) = -((high_freq_limit-F0s(F0i))./1000)./(1000./F0s(F0i))*(Cs(Ci)); %velocities in kHz/ms
	end
end

chirp_onset = zeros(1,length(vel_list));
%decode all chirp onsets
for idx = 1:length(vel_list)
	if idx == 1
		chirp_onset(idx) = 1;
	else
		chirp_onset(idx) = chirp_onset(idx-1) + chirp_npts(vel_list(idx-1)) + jitter_list(idx-1);
	end
end

% Build PSTH
if isfield(param,'stim_dur')
	%saving "params.dur" is not recommended because
	%post_process overwrites the stimulus dur with a
	%rounded duration. This segment checks if there is a
	%"stim_dur" saved and uses that as the duration (oldest
	%couple of sessions use params.dur)
	dur = param.stim_dur; % stimulus duration in seconds.
else
	dur = param.dur/1000; % stimulus duration in seconds.
end

bin_width = 1/fs_RVF;
bin_width = bin_width*1e6; % convert to microsec
stim_time_onset = stimx.times;
stim_time_offset = stim_time_onset+dur*1e6;
psth_edges = 0:bin_width:dur*1e6;
spktrains = zeros(length(psth_edges)-1,nrep);
for itime = 1:nrep
	spk_ind_this = find(cluster.t_spike >= stim_time_onset(itime) & cluster.t_spike <= stim_time_offset(itime));
	spk_time = cluster.t_spike(spk_ind_this) - stim_time_onset(itime);
	spktrains(:,itime) = histcounts(spk_time,psth_edges);
end
psth = sum(spktrains,2).';
spktrains = spktrains';

%add up spikes post chirp onsets
window_length = ((1000/min(F0s) + jitter_range(1))/1000)*fs_RVF; %window length in samples
%need to make window length a multiple of the resampling
%factor fs_RVF/new_fs
new_fs = 2000;
window_length = ceil(window_length/(fs_RVF/new_fs))*(fs_RVF/new_fs);
psth = [psth zeros(1,window_length)]; %add zeros to end of psth to make indexing later easier
spktrains = [spktrains zeros(nrep,window_length)];
rvf_psths = zeros(num_vels,window_length);

%loop through each velocity and create a psth of the time
%window equal to the maximum jitter time post chirp onset
for velidx = 1:num_vels
	this_vel_chirp_onsets = chirp_onset(vel_list==velidx);
	rep_rvf_psths = zeros(nrep,window_length);
	for i = 1:length(this_vel_chirp_onsets)
		rvf_psths(velidx,:) = rvf_psths(velidx,:) + psth(this_vel_chirp_onsets(i):this_vel_chirp_onsets(i)+window_length-1);
		rep_rvf_psths = rep_rvf_psths + spktrains(:,this_vel_chirp_onsets(i):this_vel_chirp_onsets(i)+window_length-1);
		%spktrains_by_chirp is organized like this:
		%stim_rep (full stimulus), then
		%within-SV.stim-rep aka seg_rep, then time
		spktrains_by_chirp(:,i,:) = spktrains(:,this_vel_chirp_onsets(i):this_vel_chirp_onsets(i)+window_length-1);
	end
	%resample to make peakier
	rs_rvf_psths(velidx,:) = sum(reshape(rvf_psths(velidx,:),[],length(rvf_psths(velidx,:))/(fs_RVF/new_fs)));
	rs_rep_rvf_psths = zeros(nrep,length(rs_rvf_psths));
	for itime = 1:nrep
		rs_rep_rvf_psths(itime,:) = sum(reshape(rep_rvf_psths(itime,:),[],length(rvf_psths(velidx,:))/(fs_RVF/new_fs)));
		for iseg = 1:seg_reps
			rs_spktrains_by_chirp(itime,iseg,:) = sum(reshape(spktrains_by_chirp(itime,iseg,:),[],length(rvf_psths(velidx,:))/(fs_RVF/new_fs)));
		end
	end
	%"sorta smart" method
	%find the peak of the psth, then center a 20ms window around
	%peak
	window_lim(velidx) = ceil((40/1000 + (chirp_npts(velidx)/fs_RVF))*new_fs); %this is the same windowing as in physiology
	[peak_val(velidx),peak_idx] = max(rs_rvf_psths(velidx,1:window_lim(velidx))); %find peak within valid limits (equivalent to hist_window)
	smart_window_dur = ((20/1000)*new_fs); %20 ms, in samples
	hist_window_smart(velidx,:) = [peak_idx - smart_window_dur/2 peak_idx + smart_window_dur/2]; %set window at peak +/- 10 ms
	if hist_window_smart(velidx,1) <= 0
		%if the start point is before zero
		hist_window_smart(velidx,:) = [1 1+smart_window_dur];
	elseif hist_window_smart(velidx,2) > window_lim(velidx)
		%if the end point is after the valid window
		hist_window_smart(velidx,:) = [window_lim(velidx)-smart_window_dur window_lim(velidx)];
	end
	spike_counts_smart(velidx) = sum(rs_rvf_psths(velidx,hist_window_smart(velidx,1):hist_window_smart(velidx,2)));
	rep_spike_counts(:,velidx) = sum(rs_rep_rvf_psths(:,hist_window_smart(velidx,1):hist_window_smart(velidx,2)),2);
	spike_counts_by_chirp(:,:,velidx) = sum(rs_spktrains_by_chirp(:,:,hist_window_smart(velidx,1):hist_window_smart(velidx,2)),3);
	rate_smart(velidx) = spike_counts_smart(velidx)/(seg_reps*nrep*per_reps*(smart_window_dur/new_fs)); %calculate rate: %calculate rate: num_spikes/(#SV.stim reps * #chirp_segment reps * per_reps (cycle reps in segment) * window dur)
	rep_rate(:,velidx) = rep_spike_counts(:,velidx)/(seg_reps*per_reps*(smart_window_dur/new_fs));  %remove the nrep term to get rep rate
	temp_rate_by_chirp = spike_counts_by_chirp(:,:,velidx)/(per_reps*(smart_window_dur/new_fs));
	rate_by_chirp(:,velidx) = temp_rate_by_chirp(:);
	rep_stderror(velidx) = std(rep_rate(:,velidx))/sqrt(length(rep_rate(:,velidx)));
	chirp_stderror(velidx) = std(rate_by_chirp(:,velidx))/sqrt(length(rate_by_chirp(:,velidx)));
end


headroom = 0.01;
footroom = 0.1;
graph_space = (1 - headroom - footroom)/6;
axes_height = (4/5)*graph_space;

%if we have more than 12 velocities, only plot 12 vel
%pairs' psths
if num_vels <= 12
	psths_to_plot = rs_rvf_psths;
	hist_window_to_plot = hist_window_smart;
	plot_velocities = velocities;
else
	%pick 12 vels to plot
	plot_idx = ceil(1:num_vels/12:num_vels);
	psths_to_plot = rs_rvf_psths(plot_idx,:);
	hist_window_to_plot = hist_window_smart(plot_idx,:);
	plot_velocities = velocities(plot_idx);
end
max_rate = max(max(psths_to_plot));

data.velocities = velocities;
data.rate = rate_smart;
data.rate_std = chirp_stderror;
end

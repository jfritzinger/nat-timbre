%% calculate_nat_STRF 
clear

%% Load in natural timbre file names

if ismac
	fpath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data';
else
	fpath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\';
end
tuning = readtable(fullfile(fpath, 'Tuning.xlsx')); % Load in tuning
target = 'Bassoon';

% Get all .wav files containing the target instrument name
listing = dir(fullfile(fpath, 'waveforms', ['*' target '*.wav']));
files = {listing.name};

% Extract note names and find corresponding frequencies
note_names = extractBetween(files, 'ff.', '.');
[~, index] = ismember(note_names, tuning.Note);
F0s = tuning.Frequency(index);

% Sort files and frequencies by pitch
[F0s, order] = sort(F0s);
files = files(order);

% Initialize variables for later use (if needed)
nfiles = numel(files);
wav_npts = zeros(1, nfiles);
wav_data = cell(1, nfiles);

%% Analyze natural timbre

[base, datapath, ~, ppi] = getPaths();
sheetpath = 'data/2025-manuscript/data-cleaning';
spreadsheet_name = 'PutativeTable2.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);

% Example session
putative = 'R29_TT3_P5_N02';
filename = sprintf('%s.mat', putative);
load(fullfile(datapath,'neural_data', filename)), 'data';
index = find(cellfun(@(s) strcmp(putative, s), sessions.Putative_Units));
CF = sessions.CF(index);
MTF_shape = sessions.MTF{index};

% Plot NT
params_NT = data(14, 2); % 13 is oboe, 14 is bassoon
data_NT = cell(2, 1);
if ~isempty(params_NT{1})
	data_NT = analyzeNT(params_NT{1});
end

%% 

[params_NT, fig, data] = plotNatSTRF([], params_NT, []);
%[params_NT, fig, data] = plotPhysSTRF_Simple([], params_NT, []);

%% 

function [params, fig, data] = plotPhysSTRF_Simple(cluster, params, stims)

% Process cluster
if iscell(params)
	params = params{1}; % If multiple, plots only the first instance
else
	params = params;
end
ds = params.dsid;
dsi = 1;
if isempty(stims)
	this_ds = params.stims.dsid == ds;
	cluster = params.cluster;
	stim = params.stims;
else
	if length(stims)==1
		this_ds = stims.dsid == ds(dsi);
		stim = stims;
	else
		this_ds = stims(ds).dsid == ds(dsi);
		stim = stims(ds);
	end
	params.plot_type = 'STRF';

end

% Regenerate stimuli
fs = stim.fs;
params.Fs = fs;
params.mnrep = 1;
[params] = generate_NT(params);
noises = 1:size(params.stim,1);
noise_mat = params.stim';
num_noises = length(noises);
dur = params.dur;

% Pre-event stimulus ensemble 
timerVal = tic;
win = 0.02; % 20ms time bin (20ms preceeding a spike)
T_pts = win*fs; % Number of samples for 20ms time bin 

for istim = 1:num_noises % 25 noises 
    stimi = noise_mat(:,istim).';

    bin_width = 1/fs*1e6;								% 1 sample window (convert to microsec)
    stim_time_ind = find([params.list.iwav_file] == istim);	% 10 reps, so 10 stim times per noise
    stim_time_onset = stim.times(stim_time_ind);		% stimulus time onset
    stim_time_offset = stim_time_onset+dur*1e6;			% stimulus time offset 
    psth_edges = 0:bin_width:dur*1e6;					% create bins of 1 sample each 
    spktrains = zeros(length(psth_edges)-1,length(stim_time_ind)); % create matrix [sample bins x reps]
	for itime = 1:length(stim_time_ind)					% for each rep,
		spk_ind_this = find(cluster.t_spike >= stim_time_onset(itime) & ...
			cluster.t_spike <= stim_time_offset(itime)); % get indices of all spikes that occurred in each rep
		spk_time = cluster.t_spike(spk_ind_this) - stim_time_onset(itime); % find the spike time (w.r.t. start of stimulus)
		spktrains(:,itime) = histcounts(spk_time,psth_edges'); % create psth for each rep
	end

	% PSTH for single noise 
    response = sum(spktrains,2).'; 
	% figure;
	% plot(response)
	% ylabel('Num Spikes')
	% xlabel('Time bins')

	% Calculate the spectrogram of a stimulus
	windowLength = 300; %50; % Length of the analysis window (in samples)
	noverlap = 299; %49; % Overlap size between consecutive windows (in samples)
	nfft = 960; % Length of the FFT
	[S,freq, time] = spectrogram(stimi, windowLength, noverlap, nfft, fs);
	freq(freq>16000) = [];
	S(length(freq)+1:end, :) = [];
	S = abs(S);

	% Get DC component 
	avg_stim = mean(S, 'all');

	% figure;
	% set(gcf, 'color', 'w')
	% imagesc(time, freq, S);
	% axis xy; % Flip the y-axis to match conventional spectrogram visualization
	% xlabel('Time (s)');
	% ylabel('Frequency (Hz)');
	% colorbar;
	% title(sprintf('Stimulus [%dx%d]', size(S, 1), size(S, 2)));

	% For each spike time, get spectrogram for 20ms preceding the spike
	spike_ind = find(response~=0);		% Find all spike indices 
	spike_ind(spike_ind<=T_pts) = [];	% If spike occurs before T_pts, ignore (?)
	spike_ind(spike_ind>length(time)) = [];	% If spike occurs after end of stim spectrogram, ignore (?)
	num_spikes = length(spike_ind);
	
	for ispike = 1:num_spikes
		t_spike = spike_ind(ispike);			% Time of spike
		spike = response(spike_ind(ispike));	% Number of spikes that occurred

		start_bin = t_spike-T_pts;	% Start of the time bin is 20ms before spike occurs
		end_bin = t_spike-1;			% End of time bin is one sample before spike occurs? Or time spike occurs?
		
		spectrum = S(:,start_bin:end_bin);			% Gets spectrum of stim 20ms before spike
		
		% figure;
		% t = linspace(0, win, win*fs);
		% imagesc(t, freq, spectrum)
		% axis xy; % Flip the y-axis to match conventional spectrogram visualization
		% xlabel('Time (s)');
		% ylabel('Frequency (Hz)');
		
		ST(ispike,:,:) = spectrum.*spike; 

	end
	
	% Average, making sure the responses with two spikes per bin are
	% properly averaged out 
	all_spikes = sum(response);
	avg_STRF = sum(ST,1)/all_spikes;


	avg_STRF = squeeze(avg_STRF);
	all_STRF(istim,:,:) = avg_STRF;

	% figure;
	% imagesc(avg_STRF)
	% axis xy;
	% xlabel('Samples');
	% ylabel('Freq Bins');
	% colormap(redblue)

end
elapsedTime = toc(timerVal)/60;
disp(['This took ' num2str(elapsedTime) ' minutes'])

STRF = mean(all_STRF,1);
STRF = squeeze(STRF);
STRF = flip(STRF, 2);
STRF = STRF-avg_stim;

% Scale 
t = linspace(0, win, win*fs);
f_ind = find(freq>15000|freq<300);
freq(f_ind,:) = [];
STRF(f_ind,:) = [];
flims = [300 15000];
tlims = [t(1) t(end)];
STRF_zeroed = STRF-mean(STRF);
max_STRF = max(max(abs(STRF_zeroed)));
clims_strf = max_STRF*[-1 1];

% Create figure
fig = figure;
imagesc(t, freq, STRF_zeroed, clims_strf)
%imagesc(t, f, STRF)
set(gca,'Ydir','normal', 'YLim', flims)
grid on
title('Nat Timbre STRF');
xlabel('Time (s)');
ylabel('Frequency (Hz)')
colorbar
colormap(redblue)

% Save data 
data.strf = STRF_zeroed;
data.t = t;
data.f = freq;
data.clims_strf = clims_strf;
data.tlims = tlims;
data.flims = flims;


end

%% 

function [params, fig, data] = plotNatSTRF(cluster, params, stims)

% Process cluster
if iscell(params)
	params = params{1}; % If multiple, plots only the first instance
else
	params = params;
end
ds = params.dsid;
dsi = 1;
if isempty(stims)
	this_ds = params.stims.dsid == ds;
	cluster = params.cluster;
	stim = params.stims;
else
	if length(stims)==1
		this_ds = stims.dsid == ds(dsi);
		stim = stims;
	else
		this_ds = stims(ds).dsid == ds(dsi);
		stim = stims(ds);
	end
	params.plot_type = 'STRF';

end

% Regenerate stimuli
fs = stim.fs;
params.Fs = fs;
params.mnrep = 1;
[params] = generate_NT(params);
noises = 1:size(params.stim,1);
noise_mat = params.stim';
num_noises = length(noises);
xpts = floor(params.dur*fs);
dur = params.dur;
sc = 20e-6 * power(10,params.spl/20);

% This method uses smoothed spike time
win = 0.02;
T_pts = win*fs;
h0 = zeros(1,num_noises);
h1 = zeros(num_noises,T_pts);
h2 = zeros(num_noises,T_pts,T_pts);
stim_mx = zeros(T_pts,xpts - T_pts);
window = gausswin(ceil(1e-3*fs),1.5); % 3 std

for istim = 1:num_noises
	stimi = noise_mat(:,istim).';

	bin_width = 1/fs*1e6; % convert to microsec
	stim_time_ind = find([params.list.iwav_file] == istim);
	stim_time_onset = stim.times(stim_time_ind);
	stim_time_offset = stim_time_onset+dur*1e6;
	psth_edges = 0:bin_width:dur*1e6;
	spktrains = zeros(length(psth_edges),length(stim_time_ind));
	for itime = 1:length(stim_time_ind)
		spk_ind_this = find(cluster.t_spike >= stim_time_onset(itime) & ...
			cluster.t_spike <= stim_time_offset(itime));
		spk_time = cluster.t_spike(spk_ind_this) - stim_time_onset(itime);
		spktrains(:,itime) = histc(spk_time,psth_edges);
	end
	response = sum(spktrains,2).';
	response = conv(response(T_pts + 1:end-1),window,'same');

	for T = 0:T_pts - 1
		stim_mx(T + 1,:) = stimi((T_pts + 1:end) - T);
	end
	h0(istim) = mean(response);
	h1(istim,:) = response*stim_mx';
	h2(istim,:,:) = bsxfun(@times,response - mean(response),stim_mx)*stim_mx';
end

H0 = mean(h0);
H1 = mean(h1)./sc./length(response);
H2 = reshape(mean(h2),T_pts,T_pts)./(2*sc.^2)./length(response);

if sum(H1)==0
	fig = figure;

	% Create struct that contains processed data
	data.t = [];
	data.f = [];
	data.H2ex_strf = [];
	data.H2in_strf = [];
	data.clims_strf = [];
	data.tlims = [];
	data.flims = [];
	data.strf = [];
else

	%quickplot_wk(fs, T_pts, H1, H2)
	fn = fs/2;

	tlims = [0 T_pts/fs];
	t = (0:T_pts - 1)/fs;
	f = fn*linspace(0,1,T_pts/2 + 1);

	% Process kernels and plot
	[U,S,V] = svd(H2);
	k = sign(diag(U).*diag(V)).*diag(S);	% weights for the singular vectors (see Lewis et al. 2002)
	negInd = find(k<0);						% column indices of the negatively-weighted singular vectors
	posInd = find(k>=0);					% column indices of the positively weighted singular vectors
	U = repmat(abs(k'),T_pts,1).*U;			% The weighted vectors
	U_fft = 2*abs(fft(U)/T_pts);
	U_fft = U_fft(1:T_pts/2 + 1,:);
	H2ex = U(:,posInd)*V(:,posInd)';
	H2in = U(:,negInd)*V(:,negInd)';
	H2 = H2ex + H2in;
	env = abs(hilbert(U));

	data.U = U;
	data.S = S;
	data.V = V;
	data.H0 = H0;
	data.H1 = H1;
	data.H2 = H2;
	data.T_pts = T_pts;
	data.fs = fs;

	negInd2 = negInd(1:4);
	posInd2 = posInd(1:4);

	H1_fft = 2*abs(fft(H1)/T_pts);
	H1_fft = H1_fft(1:T_pts/2 + 1);

	H2ex_strf = 2*abs(fft(H2ex)/T_pts);
	H2ex_strf = H2ex_strf(1:T_pts/2 + 1,:);
	H2ex_strf_null = median(reshape(H2ex_strf(:,t<0.0025),1,[]));
	H2ex_strf = H2ex_strf - H2ex_strf_null;

	H2in_strf = 2*abs(fft(H2in)/T_pts);
	H2in_strf = H2in_strf(1:T_pts/2 + 1,:);
	H2in_strf_null = median(reshape(H2in_strf(:,t<0.0025),1,[]));
	H2in_strf = H2in_strf - H2in_strf_null;

	clims_strf = max(max(abs([H2ex_strf,H2in_strf])))*[-1 1];

	% Get x limits
	pos_max = [db(U_fft(:,posInd2)) db(U_fft(:,negInd2))];
	yaxes_lim = -40 + max(max(db(U_fft)));
	[ind, ~] = find(pos_max > yaxes_lim);
	max_ind = max(ind);
	f_max = f(max_ind);
	flims = [0 f_max+50];

	% Plots
	fig = figure;
	tiledlayout(1, 2)
	nexttile
	axx = imagesc(t,f,H2ex_strf-H2in_strf,clims_strf);
	set(gca,'Ydir','normal','XLim',tlims,'YLim',flims)
	colormap(redblue)
	grid on
	title('Nat Timbre STRF');
	xlabel('Time (s)');
	ylabel('Frequency (Hz)')
	axx.HitTest = 'off';
	set(gca,'FontSize',7)

	nexttile
	hold on
	[~, ind] = max(db(U_fft(:,posInd2)));
	plot(f,db(U_fft(:,negInd2)),'b')
	plot(f,db(U_fft(:,posInd2)),'r')
	set(gca,'XLim',flims,'YLim',[-40 2] + max(max(db(U_fft))))
	grid on
	title(['FFT, CF = ' num2str(round(f(ind))) 'Hz']);
	xlabel('Frequency (Hz)');
	ylabel('Amplitude (dB)')
	xline(f(ind), 'k');
	set(gca,'FontSize',7)


	% Create struct that contains processed data
	data.t = t;
	data.f = f;
	data.H2ex_strf = H2ex_strf;
	data.H2in_strf = H2in_strf;
	data.clims_strf = clims_strf;
	data.tlims = tlims;
	data.flims = flims;
	data.strf = H2ex_strf-H2in_strf;
end
end
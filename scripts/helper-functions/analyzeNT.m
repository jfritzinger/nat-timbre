function data = analyzeNT(param)
% analyzeNT calculates the average rate response and the PSTH for each
%
% Inputs:
%    param - has all parameters of the stimulus and also contains the
%    cluster and stims data from the neuron.
%
% Outputs:
%    data - cell of structs that contains processed average rate data.
% 		 data{ind}.rate - average rate
%		 data{ind}.rate_std - standard deviation of rate
%		 data{ind}.pitch - pitches included in stimulus
%		 data{ind}.pitch_num - number/ordered pitches
%		 data{ind}.instrument - instruments included in stimulus
%
% Other m-files required: accumstats
%
% Depends on 'Tuning.xlsx' excel file
%
% Author: J. Fritzinger
% Created: 2024-09-16; Last revision: 2024-09-16
%
% -------------------------------------------------------------------------

ds = param.dsid;
this_ds = param.stims.dsid == ds;
cluster = param.cluster;

% Load in tuning sheet
if ismac
	tuning_dir = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/waveforms';
else
	tuning_dir = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\waveforms';
end
tuning = readtable(fullfile(tuning_dir, 'Tuning.xlsx'));

% Rate analysis
[Timbres,~,Timbrei] = unique([param.list.iwav_file].');
num_Timbres = length(Timbres);
if param.dur<1
	dur = param.dur;
else
	dur = param.dur/1000; % stimulus duration in seconds.
end
rate_size = [num_Timbres,1]; % [num_F0s,num_Fps];
spike_rates = cluster.num_spikes_delayed(this_ds)/...
	(dur - param.onsetWin/1000); % Uses onset window from post_process
if length(spike_rates)== length(Timbrei)
	[rate,rate_std] = accumstats({Timbrei},spike_rates, rate_size);

	% Sort by frequency of pitch
	index = NaN(num_Timbres, 1);
	note_names = extractBetween(param.filename, 'ff.','.wav');
	for ii = 1:num_Timbres % Find index of each note in tuning spreadsheet
		if length(note_names{ii}) > 3
			note_names{ii} = extractAfter(note_names{ii}, '.');
		end
		index(ii) = find(strcmp(note_names(ii), tuning.Note));
	end
	pitch_order = round(tuning.Frequency(index)); % Get freqs of each note
	[~, order] = sort(pitch_order); % Sort freqs
	pitch = categorical(pitch_order(order)); % Order pitches
	rate = rate(order);

	% Temporal analysis
	stim_set = find(param.stims.dsid == ds); % finds all indices for dataset
	these_spikes = ismember(cluster.abs_stim_num,stim_set); % finds all spikes that occurred in each stim set
	t_spike_rel = cluster.t_spike_rel(these_spikes); % gets timing for all spikes
	rel_id = cluster.rel_id(these_spikes); % don't know what this is 

	[vs,order2] = sort(Timbrei);
	num_stim = length(Timbrei);
	reps = zeros(num_stim,1);
	rep = zeros(num_stim,1);
	for i = 1:num_Timbres
		j = vs == i;
		reps(j) = (1:sum(j))';
	end
	rep(order2) = reps;

	for ii = 1:num_Timbres
		j = find(Timbrei == order(ii)); % finds them in proper F0 order
		x = t_spike_rel(ismember(rel_id,j));
		y = rep(rel_id(ismember(rel_id,j)));

		data.j(ii,:) = j;
	end
	data.t_spike_rel = t_spike_rel;
	data.rel_id = rel_id;
	data.rep = rep;
	

	% Calculate the predictable variance
	rate_matrix = zeros([param.nrep, num_Timbres]);
	x = reshape(spike_rates, [num_Timbres, param.nrep])';
	list_fi = reshape(Timbrei, [num_Timbres, param.nrep])';
	for irep = 1:param.nrep
		x_onerep = x(irep, :);
		fi_onerep = list_fi(irep,:);
		for j2 = 1:num_Timbres
			k = fi_onerep == j2;
			rate_matrix(irep, j2) = x_onerep(k);
		end
	end
	[V_p, ~,~] = predictableVariance(rate_matrix, Timbres);

	% Saved data
	data.rate = rate(order);
	data.rate_std = rate_std(order);
	data.F0s = pitch_order;
	data.note_names = note_names;
	data.order = order;
	data.pitch = pitch;
	data.pitch_num = pitch_order(order);
	data.V_p = V_p;
else
	data.rate = [];
	data.rate_std = [];
	data.F0s = [];
	data.note_names = [];
	data.order = [];
	data.pitch = [];
	data.pitch_num = [];
	data.V_p = NaN;
end

end
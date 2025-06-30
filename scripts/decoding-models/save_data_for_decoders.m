%% sort_neural_data.m
%
% Script that sorts all neural data for bassoon and oboe into a matrix of
%
% Author: J. Fritzinger
% Created: 2022-09-15; Last revision: 2024-10-14
%
% -------------------------------------------------------------------------
clear

%% Find natural timbre bassoon datasets, CF, and MTF

% Load in spreadsheet
addpath('/Users/jfritzinger/Projects/nat-timbre/scripts/helper-functions')
[base, datapath, savepath, ppi] = getPathsNT();
modelpath = '/Volumes/Nat-Timbre/data/manuscript';
sessions = readtable(fullfile(base, 'data-cleaning', 'Data_Table.xlsx'),...
	'PreserveVariableNames',true);

RVF_sessions = load(fullfile(base, 'RVF_PC2.mat'), "RVF_sessions");

%% Create matrices for bassoon and oboe separately

% Natural timbre datasets
NT_datasets(1,:) = cellfun(@(s) contains(s, 'R'), sessions.Oboe);
NT_datasets(2,:) = cellfun(@(s) contains(s, 'R'), sessions.Bassoon);
NT_list = find(any(NT_datasets));
num_sesh = length(NT_list);

%% Load in all data
nat_data = struct;
for ii = 1:num_sesh

	% Load in data
	putative = sessions.Putative_Units{NT_list(ii)};
	CF = sessions.CF(NT_list(ii));
	MTF_shape = sessions.MTF{NT_list(ii)};
	load(fullfile(datapath, [putative '.mat']))

	% Analyze data
	param = data(6:7, 2);
	if ~isempty(param{1})
		data_oboe = analyzeNT(param{1});
	else
		data_oboe.rate = [];
	end
	if ~isempty(param{2})
		data_bass = analyzeNT(param{2});
	else
		data_bass.rate = [];
	end

	% Get spont rate
	param_RM = data{2,2};
	if ~isempty(param_RM)
		data_RM = analyzeRM(param_RM); 
	else
		param_RM = data{2,1};
		data_RM = analyzeRM(param_RM); 
	end
	spont = data_RM.spont;

	% Get PC2 score
	has_RVF = strcmp(putative, RVF_sessions.RVF_sessions.Var1);
	if any(has_RVF)
		PC2_score = RVF_sessions.RVF_sessions.PC2_score(has_RVF);
		nat_data(ii).RVF_PC2 = PC2_score;
	else
		nat_data(ii).RVF_PC2 = NaN;
	end


	% Put all datapoints into matrix (one for BE, one BS)
	nat_data(ii).putative = putative;
	nat_data(ii).CF = CF;
	nat_data(ii).MTF = MTF_shape;

	if ~isempty(data_oboe.rate) && ~isempty(data_bass.rate)

		% Subtract spont  
		sub_rate_oboe = data_oboe.rate-spont;
		sub_rate_bass = data_bass.rate-spont;

		% Normalize rate to maximum 
		sub_rate_all = [sub_rate_oboe; sub_rate_bass];
		norm_rate = sub_rate_all / max(sub_rate_all);
		norm_rate_oboe = norm_rate(1:35);
		norm_rate_bass = norm_rate(36:end);

		% Get normalized rates per rep
		sub_rate_oboe = data_oboe.raw_rates-spont;
		sub_rate_bass = data_bass.raw_rates-spont;
		sub_rates = [sub_rate_oboe sub_rate_bass];
		max_rate = max(max(sub_rates));
		norm_rates = sub_rates / max_rate;
		norm_rates_oboe = norm_rates(:, 1:35);
		norm_rates_bass = norm_rates(:, 36:end);

	elseif  ~isempty(data_oboe.rate)
		sub_rate_oboe = data_oboe.rate-spont;
		norm_rate = sub_rate_oboe / max(sub_rate_oboe);
		norm_rate_oboe = norm_rate;

		% Get normalized rates per rep
		sub_rates = data_oboe.raw_rates-spont;
		max_rate = max(max(sub_rates));
		norm_rates_oboe = sub_rates / max_rate;

	elseif ~isempty(data_bass.rate)
		sub_rate_bass = data_bass.rate-spont;
		norm_rate = sub_rate_bass / max(sub_rate_bass);
		norm_rate_bass = norm_rate;

		% Get normalized rates per rep
		sub_rates = data_bass.raw_rates-spont;
		max_rate = max(max(sub_rates));
		norm_rates_bass = sub_rates / max_rate;
	end

	if ~isempty(data_oboe.rate)
		temporal_oboe = analyzeNT_Temporal(data_oboe, CF);
		temporal_oboe.VS_perrep(temporal_oboe.VS_perrep>=0.99) = NaN;
		temporal_oboe.VS_perrep(isnan(temporal_oboe.VS_perrep)) = 0.0001;
		sync_rep = temporal_oboe.VS_perrep' .* data_oboe.raw_rates;

		% Oboe data
		nat_data(ii).oboe_rate = data_oboe.rate;
		nat_data(ii).oboe_rate_std = data_oboe.rate_std;
		nat_data(ii).oboe_raterep = data_oboe.raw_rates;
		nat_data(ii).oboe_sync_rep = sync_rep;
		nat_data(ii).oboe_norm_rate = norm_rate_oboe; 
		nat_data(ii).oboe_norm_rep = norm_rates_oboe;
		nat_data(ii).oboe_VS = temporal_oboe.VS;
		nat_data(ii).oboe_VS_p = temporal_oboe.VS_p;
		nat_data(ii).oboe_VSrep = temporal_oboe.VS_perrep;
		nat_data(ii).oboe_VS_CF = temporal_oboe.VS_CF;
		nat_data(ii).oboe_spikerate = temporal_oboe.x;
		nat_data(ii).oboe_spikerep = temporal_oboe.y;
	end

	if ~isempty(data_bass.rate)
		temporal_bass = analyzeNT_Temporal(data_bass, CF);
		temporal_bass.VS_perrep(temporal_bass.VS_perrep>=0.99) = NaN;
		temporal_bass.VS_perrep(isnan(temporal_bass.VS_perrep)) = 0.0001;
		sync_rep = temporal_bass.VS_perrep' .* data_bass.raw_rates;

		% Bassoon data
		nat_data(ii).bass_rate = data_bass.rate;
		nat_data(ii).bass_rate_std = data_bass.rate_std;
		nat_data(ii).bass_raterep = data_bass.raw_rates;
		nat_data(ii).bass_sync_rep = sync_rep;
		nat_data(ii).bass_norm_rate = norm_rate_bass;
		nat_data(ii).bass_norm_rep = norm_rates_bass;
		nat_data(ii).bass_VS = temporal_bass.VS;
		nat_data(ii).bass_VS_p = temporal_bass.VS_p;
		nat_data(ii).bass_VSrep = temporal_bass.VS_perrep;
		nat_data(ii).bass_VS_CF = temporal_bass.VS_CF;
		nat_data(ii).bass_spikerate = temporal_bass.x;
		nat_data(ii).bass_spikerep = temporal_bass.y;
	end
	fprintf('%d%% Done!\n', round(ii/num_sesh*100))
end

%% Save dataset

save(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), "nat_data");
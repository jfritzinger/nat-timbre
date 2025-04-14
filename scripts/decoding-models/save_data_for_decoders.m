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
[base, datapath, savepath, ppi] = getPaths();
modelpath = '/Volumes/Nat-Timbre/data/manuscript';
sheetpath = 'data/2025-manuscript/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);

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
	load(fullfile(datapath, 'neural_data', [putative '.mat']))

	% Analyze data
	param = data(13:14, 2);
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

	param_RM = data{2,2};
	if ~isempty(param_RM)
		data_RM = analyzeRM(param_RM); 
	else
		param_RM = data{2,1};
		data_RM = analyzeRM(param_RM); 
	end
	spont = data_RM.spont;

	% Put all datapoints into matrix (one for BE, one BS)
	nat_data(ii).putative = putative;
	nat_data(ii).CF = CF;
	nat_data(ii).MTF = MTF_shape;

	if ~isempty(data_oboe.rate)
		temporal_oboe = analyzeNT_Temporal(data_oboe, CF);

		% Get normalized rates 
		sub_rate = data_oboe.rate-spont;
		norm_rate = sub_rate / max(sub_rate);

		% Get normalized rates per rep
		sub_rates = data_oboe.raw_rates-spont;
		max_rate = max(max(sub_rates));
		norm_rates = sub_rates / max_rate;
		% OR 
		% max_rate = max(sub_rates, [], 2);
		% norm_rates = sub_rates ./ max_rate;

		% Oboe data
		nat_data(ii).oboe_rate = data_oboe.rate;
		nat_data(ii).oboe_rate_std = data_oboe.rate_std;
		nat_data(ii).oboe_raterep = data_oboe.raw_rates;
		nat_data(ii).oboe_norm_rate = norm_rate; 
		nat_data(ii).oboe_norm_rep = norm_rates;
		nat_data(ii).oboe_VS = temporal_oboe.VS;
		nat_data(ii).oboe_VSrep = temporal_oboe.VS_rep;
		nat_data(ii).oboe_VS_CF = temporal_oboe.VS_CF;
		nat_data(ii).oboe_spikerate = temporal_oboe.x;
		nat_data(ii).oboe_spikerep = temporal_oboe.y;
	end

	if ~isempty(data_bass.rate)
		temporal_bass = analyzeNT_Temporal(data_bass, CF);

		% Get normalized rates 
		sub_rate = data_bass.rate-spont;
		norm_rate = sub_rate / max(sub_rate);

		% Get normalized rates per rep
		sub_rates = data_bass.raw_rates-spont;
		max_rate = max(max(sub_rates));
		norm_rates = sub_rates / max_rate;
		% OR 
		% max_rate = max(sub_rates, [], 2);
		% norm_rates = sub_rates ./ max_rate;

		% Bassoon data
		nat_data(ii).bass_rate = data_bass.rate;
		nat_data(ii).bass_rate_std = data_bass.rate_std;
		nat_data(ii).bass_raterep = data_bass.raw_rates;
		nat_data(ii).bass_norm_rate = norm_rate;
		nat_data(ii).bass_norm_rep = norm_rates;
		nat_data(ii).bass_VS = temporal_bass.VS;
		nat_data(ii).bass_VSrep = temporal_bass.VS_rep;
		nat_data(ii).bass_VS_CF = temporal_bass.VS_CF;
		nat_data(ii).bass_spikerate = temporal_bass.x;
		nat_data(ii).bass_spikerep = temporal_bass.y;
	end
	fprintf('%d%% Done!\n', round(ii/num_sesh*100))
end

%% Save dataset

save('Data_NT2.mat', "nat_data");





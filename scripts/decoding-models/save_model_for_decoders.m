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
%addpath('/Users/jfritzinger/Projects/nat-timbre/scripts/helper-functions')
[base, datapath, savepath, ppi] = getPathsNT();
sessions = readtable(fullfile(base, 'data-cleaning', 'Data_Table.xlsx'),...
	'PreserveVariableNames',true);

%modelpath = '/Volumes/DataFiles_JBF/Nat-Timbre/data/manuscript/SFIE_model';
modelpath = 'C:\DataFiles_JBF\Nat-Timbre\data\manuscript\SFIE_model';

%% Create matrices for bassoon and oboe separately

% Natural timbre datasets
NT_datasets(1,:) = cellfun(@(s) contains(s, 'R'), sessions.Oboe);
NT_datasets(2,:) = cellfun(@(s) contains(s, 'R'), sessions.Bassoon);
NT_list = find(any(NT_datasets));
num_sesh = length(NT_list);

%% Load in all data
nat_model = struct;
for ii = 1:num_sesh

	% Load in data
	putative = sessions.Putative_Units{NT_list(ii)};
	CF = sessions.CF(NT_list(ii));
	MTF_shape = sessions.MTF{NT_list(ii)};
	if strcmp(MTF_shape, 'BS') || strcmp(MTF_shape, 'BE')

		% Load in MTF data
		load(fullfile(modelpath, [putative '_SFIE_RM.mat']), ...
			"params_MTF", "SFIE") % Accidentally named MTF RM...
		if ~isempty(SFIE)
		nat_model(ii).MTF_rate = SFIE.rate;
		nat_model(ii).MTF_rate_std = SFIE.rate_std;
		nat_model(ii).MTF_R = SFIE.R;
		nat_model(ii).MTF_R2 = SFIE.R2;
		nat_model(ii).fms = params_MTF.all_fms;
		end

		% Load in RM data
		load(fullfile(modelpath, [putative '_SFIE_RM_actual.mat']), ...
			"params_RM", "SFIE")
		if ~isempty(SFIE)
		nat_model(ii).RM_rate = SFIE.rate;
		nat_model(ii).RM_rate_std = SFIE.rate_std;
		nat_model(ii).RM_R = SFIE.R;
		nat_model(ii).RM_R2 = SFIE.R2;
		nat_model(ii).freqs = params_RM.all_freqs;
		nat_model(ii).spls = params_RM.all_spls;
		end

		% Analyze data
		load(fullfile(modelpath, [putative '_SFIE20reps.mat']), "model_params", ...
			"params_NT", "SFIE")
		if ~isempty(SFIE{1})
			data_oboe = analyzeNTModel(SFIE{1}.SFIE_temp, params_NT{1}, MTF_shape, CF);
		else
			data_oboe.rate = [];
		end
		if ~isempty(SFIE{2})
			data_bass = analyzeNTModel(SFIE{2}.SFIE_temp, params_NT{2}, MTF_shape, CF);
		else
			data_bass.rate = [];
		end

		% Put all datapoints into matrix (one for BE, one BS)
		nat_model(ii).putative = putative;
		nat_model(ii).CF = CF;
		nat_model(ii).MTF = MTF_shape;

		if ~isempty(data_oboe.rate) % Oboe data
			nat_model(ii).oboe_rate = data_oboe.rate;
			nat_model(ii).oboe_rate_std = data_oboe.rate_std;
			nat_model(ii).oboe_raterep = data_oboe.raw_rates;
			nat_model(ii).oboe_VS = data_oboe.VS;
			nat_model(ii).oboe_VS_CF = data_oboe.VS_CF;
			nat_model(ii).oboe_PSTH = data_oboe.PSTH;
			nat_model(ii).oboe_PSTH_all =  data_oboe.PSTH_all_reps;
			nat_model(ii).oboe_reliability = data_oboe.r_splithalf;

			% Get spike train from Poisson spike generator
			x = cell(35, 1);
			y = cell(35, 1);
			for istim = 1:35
				bass_psth_one = data_oboe.PSTH_all_reps{istim};
				spike_rate = [];
				spike_rep = [];
				for irep = 1:20
					rate_fh = bass_psth_one(irep,:);
					T = 0.3;
					EventTimes = genNHPP(rate_fh,T, 1);
					spike_rate = [spike_rate EventTimes];
					num_spikes = length(EventTimes);
					spike_rep = [spike_rep; irep*ones(num_spikes, 1)];

				end
				x{istim} = spike_rate;
				y{istim} = spike_rep;
			end
			nat_model(ii).oboe_spikerate = x;
			nat_model(ii).oboe_spikerep = y;
		end

		if ~isempty(data_bass.rate) % Bassoon data
			nat_model(ii).bass_rate = data_bass.rate;
			nat_model(ii).bass_rate_std = data_bass.rate_std;
			nat_model(ii).bass_raterep = data_bass.raw_rates;
			nat_model(ii).bass_VS = data_bass.VS;
			nat_model(ii).bass_VS_CF = data_bass.VS_CF;
			nat_model(ii).bass_PSTH = data_bass.PSTH;
			nat_model(ii).bass_PSTH_all =  data_bass.PSTH_all_reps;
			nat_model(ii).bass_reliability = data_bass.r_splithalf;

			% Get spike train from Poisson spike generator
			x = cell(40, 1);
			y = cell(40, 1);
			for istim = 1:40
				bass_psth_one = data_bass.PSTH_all_reps{istim};
				spike_rate = [];
				spike_rep = [];
				for irep = 1:20
					rate_fh = bass_psth_one(irep,:);
					T = 0.3;
					EventTimes = genNHPP(rate_fh,T, 1);
					spike_rate = [spike_rate EventTimes];
					num_spikes = length(EventTimes);
					spike_rep = [spike_rep; irep*ones(num_spikes, 1)];

				end
				x{istim} = spike_rate;
				y{istim} = spike_rep;
			end
			nat_model(ii).bass_spikerate = x;
			nat_model(ii).bass_spikerep = y;
		end
	end
	fprintf('%d%% Done!\n', round(ii/num_sesh*100))
end

%% Save dataset

save(fullfile(base, 'model_comparisons', 'Model_NT3.mat'), "nat_model", '-v7.3');

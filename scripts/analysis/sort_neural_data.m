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
[base, datapath, savepath, ppi] = getPaths();
modelpath = '/Volumes/Nat-Timbre/data/manuscript';
sheetpath = 'scripts/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);

%% Create matrices for bassoon and oboe separately

for int = 1:2

	% Natural timbre datasets
	if int == 1
		NT_datasets = cellfun(@(s) contains(s, 'R'), sessions.Oboe);
	else
		NT_datasets = cellfun(@(s) contains(s, 'R'), sessions.Bassoon);
	end
	NT_list = find(NT_datasets);

	% Sort by CF
	CF_list = sessions.CF(NT_datasets);
	MTF_list = sessions.MTF(NT_datasets);
	[~, order] = sort(CF_list);
	CF_list = CF_list(order);
	MTF_list = MTF_list(order);
	NT_list = NT_list(order);
	num_data = length(NT_list);

	% Load in all data
	ind = 1;
	for ii = 1:num_data

		% Load in data
		putative = sessions.Putative_Units{NT_list(ii)};
		CF = sessions.CF(NT_list(ii));
		MTF_shape = sessions.MTF{NT_list(ii)};
		load(fullfile(datapath, [putative '.mat']))

		% Analyze data
		param = data{12+int, 2};
		data_NT = analyzeNT(param);

		% Get rid of noisy data
		if data_NT.V_p<0.5
			continue
		else
			% Put all datapoints into matrix (one for BE, one BS)
			if ~isempty(data_NT.rate)
				CF_all(ind) = CF;
				MTF_shape_all{ind} = MTF_shape;
				rate_all(ind,:) = data_NT.rate;
				rate_z_all(ind,:) = zscore(data_NT.rate);
				rate_std_all(ind,:) = data_NT.rate_std;
				ind = ind +1;
			end
		end
	end

	% Save out matrix and CFs and MTF shape
	if int == 2
		data_bassoon.CFs = CF_all;
		data_bassoon.MTF_shapes = MTF_shape_all;
		data_bassoon.rates = rate_all;
		data_bassoon.rates_z = rate_z_all;
		data_bassoon.rates_std = rate_std_all;
	else
		data_oboe.CFs = CF_all;
		data_oboe.MTF_shapes = MTF_shape_all;
		data_oboe.rates = rate_all;
		data_oboe.rates_z = rate_z_all;
		data_oboe.rates_std = rate_std_all;
		
	end
	clear CF_all MTF_shape_all rate_all rate_z_all rate_std_all
end

%% Save out matrices

save('Data_NT_Matrix.mat', "data_bassoon", "data_oboe");


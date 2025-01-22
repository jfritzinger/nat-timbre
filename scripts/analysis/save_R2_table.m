%% Instrumental Timbre Data Table 
clear 

% Load in spreadsheet
[base, datapath, savepath, ppi] = getPaths();
sheetpath = 'scripts/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);
strfpath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/strf_models';
% Initialize spreadsheet columns
modelpath = '/Volumes/Nat-Timbre/data/manuscript';
varNames = ["Putative", "CF", "MTF", "BMF"...
	"Instrument", "SFIE_R", "SFIE_R2",...
	"Energy_R", "Energy_R2", ...
	"Lat_Inh_R", "Lat_Inh_R2"];
varTypes = ["string", "double", "string", "double"...
	"string", "double", "double", ...
	"double", "double", ...
	"double", "double"];
est_num_rows = 429; % set to number larger than
num_cols = length(varNames);
table_size = [est_num_rows num_cols];
tables = table('Size',table_size,'VariableTypes',varTypes,'VariableNames',varNames);

% Find sessions for target synthetic timbre response
bin200(:,1) = cellfun(@(s) contains(s, 'R'), sessions.Bassoon);
bin200(:,2) = cellfun(@(s) contains(s, 'R'), sessions.Oboe);
isMTF = strcmp(sessions.MTF, 'BE')|strcmp(sessions.MTF, 'BS');
bin200_MTF = bin200 & isMTF;

% Add R^2 values to the spreadsheet 
has_data = bin200_MTF(:,1) | bin200_MTF(:,2);
indices = find(has_data);
num_index = length(indices);
ii = 1;
instruments = {'Oboe', 'Bassoon'};
for isesh = 1:num_index

	% Load in data
	putative = sessions.Putative_Units{indices(isesh)};
	load(fullfile(modelpath,'SFIE_model', [putative '_SFIE.mat']), 'SFIE', 'params_NT')
	load(fullfile(modelpath,'energy_model', [putative '_Energy.mat']), 'energy')
	load(fullfile(modelpath,'lat_inh_model', [putative '_Lat_Inh.mat']), 'lat_inh')
	load(fullfile(datapath, 'neural_data', [putative '.mat']))
	

	for ispl = 1:2
	
		if ispl == 1
			has_data = contains(sessions.STRF{indices(isesh)}, 'R') & ...
				contains(sessions.Oboe{indices(isesh)}, 'R');
		else
			has_data = contains(sessions.STRF{indices(isesh)}, 'R') & ...
				contains(sessions.Bassoon{indices(isesh)}, 'R');
			if strcmp(putative, 'R25_TT3_P8_N10') || strcmp(putative, 'R25_TT4_P8_N04')
				has_data = 0;
			end
		end

		if any(has_data) 
			param = data{12+ispl, 2};
			data_NT = analyzeNT(param);
			load(fullfile(strfpath, [putative '_STRF_' instruments{ispl} '.mat']))
			R = corrcoef(data_NT.rate, STRFmodel.avModel);
			STRFmodel.R = R(1,2);
		else 
			STRFmodel.R2 = NaN;
			STRFmodel.R = NaN;
		end

		if ~isempty(SFIE{ispl})
			tables.Putative{ii} = sessions.Putative_Units{indices(isesh)};
			tables.CF(ii) = sessions.CF(indices(isesh));
			tables.MTF{ii} = SFIE{ispl}.MTF_shape;
			tables.BMF(ii) = SFIE{ispl}.BMF;
			tables.Instrument{ii} = params_NT{ispl}.target;
			tables.SFIE_R(ii) = SFIE{ispl}.R;
			tables.SFIE_R2(ii) = SFIE{ispl}.R2;
			tables.Energy_R(ii) = energy{ispl}.R;
			tables.Energy_R2(ii) = energy{ispl}.R2;
			tables.Lat_Inh_R(ii) = lat_inh{ispl}.R;
			tables.Lat_Inh_R2(ii) = lat_inh{ispl}.R2;
			tables.STRF_R(ii) = STRFmodel.R;
			tables.STRF_R2(ii) = STRFmodel.R2;
			ii = ii + 1;
		end
	end

	fprintf('%s done, %d percent done\n', putative, round(isesh/num_index*100))
end

% Save table
writetable(tables,'model_r2_values_NT.xlsx')

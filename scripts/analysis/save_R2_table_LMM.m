%% Instrumental Timbre Data Table 
clear 

% Load in spreadsheet
addpath('/Users/jfritzinger/Projects/synth-timbre/scripts/helper-functions', '-end')
addpath '/Users/jfritzinger/Projects/nat-timbre/scripts/helper-functions'

[base, datapath, savepath, ppi] = getPaths();
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(datapath, 'data-cleaning', spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);
strfpath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/strf_models';
% Initialize spreadsheet columns
modelpath = '/Volumes/Nat-Timbre/data/manuscript';
varNames = ["Putative", "CF", "MTF", "BMF"...
	"Instrument", "F0", "Model", "R"];
varTypes = ["string", "double", "string", "double"...
	"string", "string", "string", "double"];
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
model_names = {'SFIE', 'BroadInh', 'Energy'};
for isesh = 1:num_index

	% Load in data
	putative = sessions.Putative_Units{indices(isesh)};
	load(fullfile(modelpath,'SFIE_model', [putative '_SFIE.mat']), 'SFIE', 'params_NT')
	load(fullfile(modelpath,'energy_model', [putative '_Energy.mat']), 'energy')
	load(fullfile(modelpath,'lat_inh_model', [putative '_Lat_Inh.mat']), 'lat_inh')
	load(fullfile(datapath, 'neural_data', [putative '.mat']))

	for iinst = 2 % Only bassoon

		% Add voice and high F0 R^2 seperately
		param = data{12+iinst, 2};
		for imodel = 1:3
			for iF0 = 1:2
				if ~isempty(param) && ~isempty(SFIE{iinst})
					data_NT = analyzeNT(param);

					if iF0 == 1
						index = data_NT.pitch_num<260 & data_NT.pitch_num>80;
						F0 = 'Voice';
					else
						index = data_NT.pitch_num>260;
						F0 = 'High';
					end

					if imodel == 1 % SFIE
						R = corrcoef(data_NT.rate(index), SFIE{iinst}.rate(index));
					elseif imodel == 2 % Broad inhibition
						R = corrcoef(data_NT.rate(index), lat_inh{iinst}.rate(index));
					else % Energy
						R = corrcoef(data_NT.rate(index), energy{iinst}.rate(index));
					end
				else
					R = [NaN NaN];
				end

				if ~isempty(SFIE{iinst})
					tables.Putative{ii} = sessions.Putative_Units{indices(isesh)};
					tables.CF(ii) = sessions.CF(indices(isesh));
					tables.MTF{ii} = SFIE{iinst}.MTF_shape;
					tables.BMF(ii) = SFIE{iinst}.BMF;
					tables.Instrument{ii} = params_NT{iinst}.target;
					tables.F0{ii} = F0;
					tables.Model{ii} = model_names{imodel};
					tables.R(ii) = R(1,2);
					ii = ii + 1;
				end
			end
		end
	end
	fprintf('%s done, %d percent done\n', putative, round(isesh/num_index*100))
end

% Save table
writetable(tables,'model_r2_LMM.xlsx')

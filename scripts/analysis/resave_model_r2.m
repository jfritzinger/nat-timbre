%% save_model_predictions_NT
%
% Script that runs either SFIE single-cell, energy, or SFIE population
% models for each neuron with a response to synthetic timbre
%
% Author: J. Fritzinger
% Created: 2022-09-15; Last revision: 2024-09-25
%
% -------------------------------------------------------------------------
clear

%%

%model_type = 'SFIE';
%model_type = 'Energy';
model_type = 'Lat_Inh';

% Load in spreadsheet
[base, datapath, savepath, ppi] = getPaths();
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(datapath, 'data-cleaning', spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);

% Find sessions for target synthetic timbre response
timbre(:,1) = cellfun(@(s) contains(s, 'R'), sessions.Bassoon);
timbre(:,2) = cellfun(@(s) contains(s, 'R'), sessions.Oboe);

%%

has_data = timbre(:,1) | timbre(:,2);
indices = find(has_data);
num_index = length(indices);
for isesh = 1:num_index
	timerVal = tic;

	% Load in session
	putative = sessions.Putative_Units{indices(isesh)};
	CF = sessions.CF(indices(isesh));
	load(fullfile(datapath, 'neural_data', [putative '.mat']))
	MTF_shape = sessions.MTF{indices(isesh)};
	params_NT = data(13:14, 2); %%%%%%%%%%%%%%%

	% Load model
	filename = [putative '_' model_type '.mat'];
	filepath = 'C:\DataFiles_JBF\Nat-Timbre\data\manuscript';
	%filepath = '/Volumes/Nat-Timbre/data/manuscript/';
	if strcmp(model_type, 'Energy')
		load(fullfile(filepath, 'energy_model', filename), 'params_NT', 'energy')
	elseif strcmp(model_type, 'SFIE')
		load(fullfile(filepath, 'SFIE_model', filename), 'params_NT', 'AN', 'SFIE', 'model_params')
	elseif strcmp(model_type, 'Lat_Inh')
		load(fullfile(filepath, 'lat_inh_model', filename), 'params_NT', 'AN_lat_inh', 'lat_inh', 'model_params')
	end

	for iNT = 1:2

		% Analysis
		if isempty(params_NT{iNT})
			% Did not record natural timbre response for this instrument
		else
			data_NT = analyzeNT(params_NT{iNT});
			if isempty(data_NT.rate)

			else
				if strcmp(model_type, 'SFIE')
					if strcmp(MTF_shape, 'BS') || strcmp(MTF_shape, 'BE')
						R = corrcoef(data_NT.rate,SFIE{iNT}.rate); % Correlation
						SFIE{iNT}.R = R(1, 2);
						SFIE{iNT}.R2 = R(1, 2).^2;
					end
					R = corrcoef(data_NT.rate,AN{iNT}.rate);
					AN{iNT}.R = R(1, 2);
					AN{iNT}.R2 = R(1, 2).^2;

				elseif strcmp(model_type, 'Energy') % Energy model
					R_int = corrcoef(data_NT.rate,energy{iNT}.rate);
					energy{iNT}.R = R_int(1,2);
					energy{iNT}.R2 =  R_int(1, 2).^2;

				elseif strcmp(model_type, 'Lat_Inh')

					if strcmp(MTF_shape, 'BS') || strcmp(MTF_shape, 'BE')
						R = corrcoef(data_NT.rate,lat_inh{iNT}.rate); % Correlation
						lat_inh{iNT}.R = R(1, 2);
						lat_inh{iNT}.R2 = R(1, 2).^2;
					end
					R = corrcoef(data_NT.rate,AN_lat_inh{iNT}.rate);
					AN_lat_inh{iNT}.R = R(1, 2);
					AN_lat_inh{iNT}.R2 = R(1, 2).^2;
				end	
			end
		end
	end

	% Save model
	filename = [putative '_' model_type '.mat'];
	if strcmp(model_type, 'Energy')
		save(fullfile(filepath, 'energy_model', filename), 'params_NT', 'energy')
	elseif strcmp(model_type, 'SFIE')
		save(fullfile(filepath, 'SFIE_model', filename), 'params_NT', 'AN', 'SFIE', 'model_params')
	elseif strcmp(model_type, 'Lat_Inh')
		save(fullfile(filepath, 'lat_inh_model', filename), 'params_NT', 'AN_lat_inh', 'lat_inh', 'model_params')
	end

	elapsedTime = toc(timerVal);
	disp([putative ' Model took ' num2str(elapsedTime) ' seconds'])
end

%% Create putative neuron spreadsheet 
% J. Fritzinger
%
clear 

%% Load in spreadsheet of neurons that have putative units 

spreadsheet_name = 'NTSessions_2025.xlsx';
sessions = readtable(spreadsheet_name, 'PreserveVariableNames',true);
sessions(strcmp(sessions.Putative_Units,""),:) = [];

% Reads in spreadsheet & load data
% if ismac
% 	filename = '/Volumes/Rabbit_data/session_table.xlsx';
% else
% 	filename = '\\nsc-lcarney-g1\Rabbit_data\session_table.xlsx';
% end
% cellData = readtable(filename, 'PreserveVariableNames',true);

%% Find all putative units 

putative = sessions.Putative_Units;
[labels_putative, I, C] = unique(putative);
num_putative_neurons = length(labels_putative);

%% Set up spreadsheet 

% Initialize spreadsheet columns - work in progress
varNames = ["Putative_Units", "CF", "MTF", ...
	"BMF", "WMF", "MTF_con", "BMF_con", "WMF_con"];
varTypes = ["string", "double", "string", ...
	"double", "double", "string", "double", "double"];

stimNames = ["char_spl", "char_ITD", "char_ILD", "type=RM", "type=RM_con", ...
	"typMTFN","typMTFN_con","STRF","STRF_con","SCHR" ...
	"43dB", "63dB", "73dB", "83dB", "43dB_100", "63dB_100", "73dB_100", "83dB_100",...
	"43dB_con", "63dB_con", "73dB_con", "83dB_con","Oboe", "Bassoon", "Other"];

stimTypes = repmat("string", 1, length(stimNames));
est_num_rows = 75; % set to number larger than
num_cols = length([varNames stimNames]);
table_size = [est_num_rows num_cols];
putative_table = table('Size',table_size,'VariableTypes',[varTypes stimTypes],'VariableNames',[varNames stimNames]);

for ind = 1:num_putative_neurons 

	num_repeats = sum(C==ind);
	putative_table.Putative_Units(ind) = labels_putative(ind);
	if num_repeats == 1

		table_index = find(C==ind);
		rabbit = sessions.Rabbit(table_index);
		session = sessions.Session(table_index);
		tetrode = sessions.TT(table_index);
		neuron = sessions.N(table_index);

		putative_table(ind, 2:8) = sessions(table_index, 10:16); % change

		% Basic stimuli: 9-18 are basic stimuli
		for ids = 18:27
			if sessions{table_index, ids} == 1
				name = sprintf('R0%dS%3g_TT%d_N%d', rabbit, session, tetrode, neuron);
				name = regexprep(name, ' ', '0');
				putative_table(ind, ids-9) = {name};
			end
		end

		% Natural timbre: 30-32 are natural timbre
		if strcmp(sessions.Include_NT{ind}, "Y")
			for ids = 40:42
				if sessions{table_index, ids} == 1
					name = sprintf('R0%dS%3g_TT%d_N%d', rabbit, session, tetrode, neuron);
					name = regexprep(name, ' ', '0');
					putative_table(ind, ids-9) = {name};
				end
			end
		end
	end
end

%% How many to do:

missing = sum(putative_table.CF==0);

%% Save Putative Table 

writetable(putative_table, 'Data_Table_2025.xlsx')
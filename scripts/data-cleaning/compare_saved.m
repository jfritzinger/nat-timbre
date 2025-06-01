clear  

%% get putative units 

spreadsheet_name = 'Data_Table.xlsx';
sessions = readtable(spreadsheet_name, 'PreserveVariableNames',true);
num_data = size(sessions, 1);

table_putative = sessions.Putative_Units;

%% get all saved files 

if ismac
	path = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/neural_data';
else
	path = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\neural_data';
end

data = dir(path);
data_putative = {data(3:end).name}';


%% compare

for idata = 1:299
	putative2 = data_putative{idata};

	for itable = 1:299

		putative1 = table_putative{itable};
		match(itable) = contains(putative2,putative1);

	end

	if sum(match) == 0
		error('Check!')
	elseif sum(match) == 1
		disp('Good')
	else
		error('Something wrong!')
	end
end


%% Wrong: 
%% 

%% Load in excel file of sessions 

% Reads in spreadsheet 
filename = 'NTSessions_All.xlsx';
sessions = readtable(filename, 'PreserveVariableNames',true);

%% Cut down spreadsheet to just 'Y's

no_SC = strcmp(sessions.Include_SC, 'N');
no_NT = strcmp(sessions.Include_NT, 'N');
wrongF0 = strcmp(sessions.Error, 'Wrong F0');
no_both = no_SC & no_NT;
num_not = sum(no_both);

sessions_new = sessions(~no_both,:);


%% Save spreadsheet 

writetable(sessions_new,'NTSessions_2025.xlsx')
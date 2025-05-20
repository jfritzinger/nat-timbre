%% Fig0_DataToReanalyze
% J. Fritzinger
%
%
clear 

%% Load in excel file of sessions 

% Reads in spreadsheet 
if ismac
	mount('nsc-lcarney-g1')
    filename = '/Volumes/Rabbit_data/session_table.xlsx';
else
	filename = '\\nsc-lcarney-g2\Rabbit_data\session_table.xlsx';
end
sessions = readtable(filename, 'PreserveVariableNames',true);
sessions(sessions.Rabbit == 0,:) = [];

%% Do analysis only for data after these session numbers 

interest_27 = sessions.Rabbit==27 & sessions.Session>421;
interest_29 = sessions.Rabbit==29 & sessions.Session>358;

%% Step 1: Find all sessions that haven't been analyzed at all 
% Print all sessions with errors 
% 
% errorind = ~strcmp(sessions.Error_Column,'') & (interest_27 | interest_29);
% 
% sessions_to_reanalyze = sessions.Session(errorind);
% rabbits_to_reanalyze = sessions.Rabbit(errorind);
% 
% for ii = 1:length(sessions_to_reanalyze)
% 	name{ii} = sprintf('R0%gS%3g', rabbits_to_reanalyze(ii), sessions_to_reanalyze(ii));
% end
% name = unique(name)';
% list = table(name);
% 
% outputfilename = 'ErroredSessions.xlsx';
% writetable(list,outputfilename)
% disp('Successfully created spreadsheet!')

%% Step 2: Print list of all sessions with timbre

stim2ind = sessions.Natural_Timbre == 1;

index = stim2ind & (interest_27 | interest_29);

sessions_to_reanalyze = sessions.Session(index);
rabbits_to_reanalyze = sessions.Rabbit(index);

for ii = 1:length(sessions_to_reanalyze)
	name{ii} = sprintf('R0%gS%3g', rabbits_to_reanalyze(ii), sessions_to_reanalyze(ii));
end
name = unique(name);
name = name';

% Save list 
list = table(name);
outputfilename = 'DataToReanalyze.xlsx';
writetable(list,outputfilename)
disp('Successfully created spreadsheet!')

%% Step 3: Get sessions that need to have CF and tetrode depth added 

rating_ind = strcmp(sessions.Rating, 'Excellent');
rating2_ind = strcmp(sessions.Rating, 'Good');
stim2ind = sessions.Natural_Timbre == 1;

index = (rating_ind|rating2_ind) & (stim2ind) & (interest_27 | interest_29);
CF_missing = isnan(sessions.CF);
Depth_missing = isnan(sessions.Tetrode_Depth);
metaind = index & (CF_missing | Depth_missing); 

sessions_to_reanalyze = sessions.Session(metaind);
rabbits_to_reanalyze = sessions.Rabbit(metaind);
for ii = 1:length(sessions_to_reanalyze)
	name{ii} = sprintf('R0%gS%3g', rabbits_to_reanalyze(ii), sessions_to_reanalyze(ii));
end
CFs = sessions.CF(metaind);
depths = sessions.Tetrode_Depth(metaind);


list = table(name', CFs, depths);

outputfilename = 'AddMetadata.xlsx';
writetable(list,outputfilename)
disp('Successfully created spreadsheet!')

%% Neuron Counts 

rating_ind = strcmp(sessions.Rating, 'Excellent');
rating2_ind = strcmp(sessions.Rating, 'Good');
stim2ind = sessions.Natural_Timbre == 1;
notdriven = sessions.CF==0;

index = (rating_ind|rating2_ind) & (stim2ind) & (interest_27 | interest_29) & ~notdriven;
sessions_nattimbre = sessions(index,:);

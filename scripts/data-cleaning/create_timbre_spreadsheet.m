%% New Spreadsheet with all Spectral Centroid data 
% J. Fritzinger, updated 1/11/24
clear

%% Load in stimuli_sessions 

% Reads in spreadsheet & load data
if ismac
	filename = '/Volumes/Rabbit_data/session_table.xlsx';
else
	filename = '\\nsc-lcarney-g1\Rabbit_data\session_table.xlsx';
end
sessions = readtable(filename, 'PreserveVariableNames',true);

%% Find sessions of interest 


rating = {'Good', 'Excellent'};
interest_rating = ismember(sessions.Rating,rating);
driven = sessions.CF~=0;

interest_27 = sessions.Rabbit==27 & sessions.Session>421;
interest_29 = sessions.Rabbit==29 & sessions.Session>358;

stimulus2 = 'Natural_Timbre';
names = sessions.Properties.VariableNames;
ind = find(strcmp(names,stimulus2));
stim_columns2 = table2array(sessions(:,ind));
interest_stim2 = sum(stim_columns2,2);


interest = (interest_27 | interest_29) & ...
	interest_rating & interest_stim2 & driven;
index = find(interest);
num_interest = length(index);


%% Create empty table 

varNames = ["Rabbit", "Session","TT","N", "Putative_Units", ...
	"Include_SC", "Include_NT", "Pass", "Depth", "CF", "MTF", ...
	"BMF", "WMF", "MTF_con", "BMF_con", "WMF_con", "Error"];
varTypes = ["double", "double","double","double", "string", ...
	"string","string", "double", "double", ...
	"double", "string", "double", "double",...
	"string", "double", "double", "string"];
stimNames = ["char_spl", "char_ITD", "char_ILD", "type_RM", "type_RM_con", ...
	"typMTFN","typMTFN_con","STRF","STRF_con","SCHR","RVF", ...
	"Oboe", "Bassoon", "Other"];

stimTypes = repmat("double", 1, length(stimNames));
est_num_rows = 1000; % set to number larger than
num_cols = length([varNames stimNames]);
table_size = [est_num_rows num_cols];
listData = table('Size',table_size,'VariableTypes',[varTypes stimTypes],'VariableNames',[varNames stimNames]);


%% Fill out table (from post_process)

[userid, base_dir, ~, report_path, data_path] = findPaths();
for iclus = 1:num_interest

	% Fill out table from sessions 
	session = sprintf('R%03dS%03d', sessions.Rabbit(index(iclus)), sessions.Session(index(iclus)));
	
	listData.Rabbit(iclus) = sessions.Rabbit(index(iclus));
	listData.Session(iclus) = sessions.Session(index(iclus));
	listData.TT(iclus) = sessions.Tetrode(index(iclus));
	listData.N(iclus) = sessions.Neuron(index(iclus));
	listData.Pass(iclus) = sessions.Pass(index(iclus));
	listData.Depth(iclus) = sessions.Tetrode_Depth(index(iclus));

	listData.CF(iclus) = sessions.CF(index(iclus));
	listData.MTF(iclus) = sessions.MTF(index(iclus));
	listData.BMF(iclus) = sessions.BMF(index(iclus));
	listData.WMF(iclus) = sessions.WMF(index(iclus));
	listData.STRF(iclus) = sessions.STRF(index(iclus));
	listData.char_spl(iclus) = sessions.char_spl(index(iclus));
	listData.char_ITD(iclus) = sessions.char_ITD(index(iclus));
	listData.char_ILD(iclus) = sessions.char_ILD(index(iclus));
	listData.type_RM(iclus) = sessions.type_RM(index(iclus));
	listData.typMTFN(iclus) = sessions.typMTFN(index(iclus));
	listData.RVF(iclus) = sessions.RVF(index(iclus));
	listData.SCHR(iclus) = sessions.SCHR(index(iclus));

	% Get the path to each session
	rab_num = num2str(sessions.Rabbit(index(iclus)));
	CF = sessions.CF(index(iclus));
	if ismac
		session_dir_name = fullfile(base_dir, ['R0' num2str(rab_num)]);
	else
		session_dir_name = base_dir{contains(base_dir, rab_num)};
	end
	session_dir = fullfile(session_dir_name, session);

	% Load in session, or add message saying it could not be loaded
	try
	[clusters, params, ~] = loadPhysiologySession(session_dir, session, userid);
	cluster = clusters([clusters.tetrode] == sessions.Tetrode(index(iclus))); % Select tetrode
	cluster = cluster([cluster.neuron] == sessions.Neuron(index(iclus))); % Select neuron
	catch
		continue
	end

	 % Contra 
	 RM_con = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'type=RM') &&...
		 p.binmode==1, params);
	listData.type_RM_con(iclus) = any(RM_con);
	MTFN_con = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'typMTFN') &&...
		 p.binmode==1, params);
	listData.typMTFN_con(iclus) = any(MTFN_con);
	STRF_con = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'STRF') &&...
		 p.binmode==1, params);
	listData.STRF_con(iclus) = any(STRF_con);

	% Get data for rest of table 
	[oboe, bassoon, other] = finddata(params);
	listData.("Oboe")(iclus) = any(oboe);
	listData.("Bassoon")(iclus) = any(bassoon);
	listData.("Other")(iclus) = any(other);

end

%% Save spreadsheet

writetable(listData,'NTSessions_All2.xlsx')


%% Functions 

function [oboe, bassoon, other] = finddata(population)

	% Find oboe 
	oboe = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'Natural_Timbre') && ...
		(isfield(p, 'target') && strcmp(p.target, 'Oboe')||...
		contains(p.list(1).wav_file, 'Oboe')), population);

	% Find bassoon
	bassoon = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'Natural_Timbre') && ...
		(isfield(p, 'target') && strcmp(p.target, 'Bassoon')||...
		contains(p.list(1).wav_file, 'Bassoon')), population);

	% Find other
	other = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'Natural_Timbre'), population);

end
%% energy_decoding
clear


%% Find natural timbre bassoon datasets, CF, and MTF

% Load in spreadsheet
addpath('/Users/jfritzinger/Projects/nat-timbre/scripts/helper-functions')
[base, datapath, savepath, ppi] = getPaths();
sheetpath = 'data/2025-manuscript/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name),...
	'PreserveVariableNames',true);

%% Create matrices for bassoon and oboe separately

% Natural timbre datasets
NT_datasets(1,:) = cellfun(@(s) contains(s, 'R'), sessions.Oboe);
NT_datasets(2,:) = cellfun(@(s) contains(s, 'R'), sessions.Bassoon);
NT_list = find(any(NT_datasets));
num_sesh = length(NT_list);
[base, ~, ~, ~] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Energy_NT.mat'), 'nat_data')

%% Load in all data
% % Load in energy
% if ismac
% 	modelpath = '/Volumes/Nat-Timbre/data/manuscript';
% else
% 	modelpath = 'C:\DataFiles_JBF\Nat-Timbre\data\manuscript';
% end
% nat_data = struct;
% for ii = 1:num_sesh
% 
% 	% Load in data
% 	putative = sessions.Putative_Units{NT_list(ii)};
% 	CF = sessions.CF(NT_list(ii));
% 	load(fullfile(modelpath,'energy_model', [putative '_Energy.mat']), 'energy')
% 
% 	% Put all datapoints into matrix (one for BE, one BS)
% 	nat_data(ii).putative = putative;
% 	nat_data(ii).CF = CF;
% 
% 	if ~isempty(energy{1})
% 		nat_data(ii).oboe_rate = energy{1}.rate;
% 	end
% 
% 	if ~isempty(energy{2})
% 		nat_data(ii).bass_rate = energy{2}.rate;
% 	end
% 	fprintf('%d%% Done!\n', round(ii/num_sesh*100))
% end
% 
% % Save output 
% save('Energy_NT.mat', 'nat_data')

%% Reshape matrix

%target = 'Bassoon';
target = 'Oboe';

if ismac
	fpath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data';
else
	fpath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\';
end
tuning = readtable(fullfile(fpath, 'Tuning.xlsx')); % Load in tuning

% Get bassoon stimulus
listing = dir(fullfile(fpath, 'waveforms', ['*' target '*.wav']));
files = {listing.name};
note_names = extractBetween(files, 'ff.', '.');
[~, index] = ismember(note_names, tuning.Note);
F0s1 = tuning.Frequency(index);
[F0s, order] = sort(F0s1);

% Find all rows with bassoon in them
%sesh = find(~cellfun(@isempty, {nat_data.bass_rate}));
sesh = find(~cellfun(@isempty, {nat_data.oboe_rate}));
num_data = numel(sesh);

data_mat = NaN(length(F0s)*20, num_data);
for ii = 1:num_data
	X1 = nat_data(sesh(ii)).oboe_rate';
	X2 = reshape(repmat(X1', 1, 20)', 1, []);
	data_mat(:,ii) = X2;
end

% Create array of correct responses
response = reshape(repmat(F0s, 1, 20)', 1, []);

% Create table for model
T = array2table(data_mat);
T.Response = response';

% % Use only first 20 F0s 
% T_new = T(1:400, :);
% 
% % Use only odd F0d
% idx = [];
% for k = 0:19
%     idx = [idx, (1 + 40*k):(20 + 40*k)];
% end
% T_new2 = T(idx,:);

%% Model - same type of model as data 

[trainedClassifier, validationAccuracy, validationPredictions] = ...
	trainClassifierPopRateF0(T, F0s);

%% Plot confusion matrix 

figure

nexttile
C = confusionmat(T.Response, validationPredictions);
confusionchart(C)
title(['Energy, ' target ', Accuracy = ' num2str(round(validationAccuracy*100)) '%'])

%% model_pop_VS_F0
clear

%% Load in data

[base, ~, ~, ~] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')

%% Get correct output of model
target = 'Bassoon';
%target = 'Oboe';

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

% Find all rows with bassoon and oboe
if strcmp(target, 'Bassoon')
	sesh_all = find(~cellfun(@isempty, {nat_data.bass_rate}));
else
	sesh_all = find(~cellfun(@isempty, {nat_data.oboe_rate}));
end
num_data_all = numel(sesh_all);
CFs_all = [nat_data(sesh_all).CF];
putative = {nat_data(sesh_all).putative};
MTFs = {nat_data(sesh_all).MTF};

%% Run model

CF_groups = [0, 14000; 0, 2000; 2000, 4000; 4000, 14000];
%accuracies = NaN(7, 500);
%totalnum = 88; % bassoon
totalnum = 30;
for iCF = 1:7
	timerVal = tic;
	parfor irep = 1:1000

		% Get random assortment of 60 units from each (for all, gets 20 per
		% group)
		ind_low = find(CFs_all > CF_groups(2, 1) & CFs_all < CF_groups(2, 2));
		ind_med = find(CFs_all > CF_groups(3, 1) & CFs_all < CF_groups(3, 2));
		ind_high = find(CFs_all > CF_groups(4, 1) & CFs_all < CF_groups(4, 2));
		if iCF == 1
			rand_ind = [randsample(ind_low, int8(totalnum/3), false) ...
				randsample(ind_med, int8(totalnum/3), false)...
				randsample(ind_high, int8(totalnum/3), false)];
		elseif ismember(iCF, [2, 3, 4])
			ind = CFs_all > CF_groups(iCF, 1) & CFs_all < CF_groups(iCF, 2);
			rand_ind = randsample(find(ind), totalnum, false);
		elseif iCF == 5
			rand_ind = [randsample(ind_low, totalnum/2, false) ...
				randsample(ind_med, totalnum/2, false)];
		elseif iCF == 6
			rand_ind = [randsample(ind_low, totalnum/2, false) ...
				randsample(ind_high, totalnum/2, false)];
		elseif iCF == 7
			rand_ind = [randsample(ind_med, totalnum/2, false) ...
				randsample(ind_high, totalnum/2, false)];
		end

		% Subset!
		CFs = CFs_all(rand_ind);
		sesh = sesh_all(rand_ind);
		num_data = length(sesh);

		% Model including all F0s
		% Find all rows with bassoon in them
		if strcmp(target, 'Bassoon') % Bassoon
			data_mat = NaN(length(F0s)*20, num_data);
			for ii = 1:num_data
				X1 = nat_data(sesh(ii)).bass_raterep';
				X2 = reshape(X1, [], 1);
				data_mat(:,ii) = X2;
			end
		else % Oboe
			data_mat = NaN(length(F0s)*20, num_data);
			for ii = 1:num_data
				X1 = nat_data(sesh(ii)).oboe_raterep';
				X2 = reshape(X1, [], 1);
				data_mat(:,ii) = X2;
			end
		end
		num_data = numel(sesh);

		% Create table for model
		response = reshape(repmat(F0s, 1, 20)', 1, []);
		T = array2table(data_mat);
		T.Response = response';

		% Run model
		if strcmp(target, 'Bassoon')
			[trainedClassifier1, validationAccuracy1, validationPredictions1] = ...
				trainClassifierPopRateF0Bass(T);
		else
			[trainedClassifier1, validationAccuracy1, validationPredictions1] = ...
				trainClassifierPopRateF0Oboe(T);
		end
		C = confusionmat(response, validationPredictions1);

		accuracies(iCF, irep) = sum(diag(C))/sum(C, 'all');
		%fprintf('%d/%d, %0.2f%% done!\n', irep, 500, irep/500*100)
	end

	% Print out progress
	timer = toc(timerVal);
	fprintf('Models took %0.2g minutes\n', timer/60)
	fprintf('%d/%d, %0.2f%% done!\n', iCF, 7, iCF/7*100)
end

%% 
CF_names = {'Low', 'Medium', 'High', 'Low+Med', 'Low+High', 'Med+High','All'};

figure
boxchart(accuracies')
xticklabels(CF_names)
hold on
max_accur = max(accuracies, [], 2);
plot(1:7, max_accur)


%% Plot outputs
% 
% save(fullfile(base, 'model_comparisons', 'Pop_Rate_F0_Subset_CF.mat'), ...
% 	"accuracies")


%% model_pop_timbre_subsets

%% model_pop_VS_F0
clear

%% Load in data

[base, ~, ~, ~] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')

%% Get correct output of model

tuning = readtable(fullfile(base, 'Tuning.xlsx')); % Load in tuning
target = 'Bassoon';
listing = dir(fullfile(base, 'waveforms', ['*' target '*.wav']));
files = {listing.name};
note_names = extractBetween(files, 'ff.', '.');
[~, index] = ismember(note_names, tuning.Note);
F0s1 = tuning.Frequency(index);
[F0s, order] = sort(F0s1);

%% Get data into proper matrix
type = 'rate';

figure
tiledlayout(4, 2)

% Find all rows with bassoon and oboe
has_bass = ~cellfun(@isempty, {nat_data.bass_rate});
has_oboe = ~cellfun(@isempty, {nat_data.oboe_rate});
sesh_all = find(has_bass & has_oboe);
num_data_all = numel(sesh_all);
ind_b = 25:40;
ind_o = [1 3:17];
CFs_all = [nat_data(sesh_all).CF];
putative = {nat_data(sesh_all).putative};
MTFs = {nat_data(sesh_all).MTF};

% ind_all(1,:) = strcmp('BE', MTFs);
% ind_all(2,:) = strcmp('BS', MTFs);
% 
% CFs_all = CFs_all(ind_all(1,:));
% putative = putative(ind_all(1,:));
% MTFs = MTFs(ind_all(1,:));

% Get a subset of the data based on CF grouping
CF_groups = [0, 14000; 0, 2000; 2000, 4000; 4000, 14000];
for iCF = 1:7

	for irep = 1:500
		
		% Get random assortment of 60 units from each (for all, gets 20 per
		% group)
		if iCF == 1
			ind_low = find(CFs_all > CF_groups(2, 1) & CFs_all < CF_groups(2, 2));
			ind_med = find(CFs_all > CF_groups(3, 1) & CFs_all < CF_groups(3, 2));
			ind_high = find(CFs_all > CF_groups(4, 1) & CFs_all < CF_groups(4, 2));
			rand_ind = [randsample(ind_low, 20, false) ...
				randsample(ind_med, 20, false)...
				randsample(ind_high, 20, false)];
		elseif ismember(iCF, [2, 3, 4])
			ind = CFs_all > CF_groups(iCF, 1) & CFs_all < CF_groups(iCF, 2);
			rand_ind = randsample(find(ind), 60, false);
		elseif iCF == 5
			rand_ind = [randsample(ind_low, 30, false) ...
				randsample(ind_med, 30, false)];
		elseif iCF == 6
			rand_ind = [randsample(ind_low, 30, false) ...
				randsample(ind_high, 30, false)];
		elseif iCF == 7
			rand_ind = [randsample(ind_med, 30, false) ...
				randsample(ind_high, 30, false)];
		end


		% Get a subset of the data based on MTF type
		% isBE = strcmp('BE', MTFs);
		% isBS = strcmp('BS', MTFs);
		% isH = contains(MTFs,'H');
		% isF = strcmp('F', MTFs);
		% ind = isF;

		% Subset!
		CFs = CFs_all(rand_ind);
		sesh = sesh_all(rand_ind);
		num_data = numel(sesh);

		%% Model including all F0s

		% Create array of correct responses
		response = [ones(1, 20) repmat(2, 1, 20)]';

		data_mat = NaN(2*20, num_data);
		data_mat2 = NaN(2*20*16, num_data);
		for target = 1:16
			for ii = 1:num_data

				% Arrange data for SVM
				X1 = nat_data(sesh(ii)).bass_raterep(:,ind_b(target));
				X2 = nat_data(sesh(ii)).oboe_raterep(:,ind_o(target));
				X = [X1; X2];
				data_mat(:,ii) = X;
			end
			idx = (1:40) + 40*(target-1);
			data_mat2(idx, :) = data_mat;
		end

		% Put data into table
		T = array2table(data_mat2);
		T.Instrument = repmat([ones(20,1); ones(20, 1)*2], 16, 1);

		%% Run model with kfold validation

		[trainedClassifier, accuracy, predictions] = trainClassifierPopRateTimbre(T);
		pop_rate_timbre.trainedClassifier = trainedClassifier;
		pop_rate_timbre.accuracy = accuracy;
		pop_rate_timbre.predictions = predictions;
		pop_rate_timbre.response = T.Instrument;
		pop_rate_timbre.T = T;
		pop_rate_timbre.CFs = CFs;
		pop_rate_timbre.putative = putative;
		pop_rate_timbre.sesh = sesh;
		pop_rate_timbre.MTF = {nat_data(sesh).MTF};
		pop_rate_timbre.oboe_rate = [nat_data(sesh).oboe_rate];
		pop_rate_timbre.oboe_rate_std = [nat_data(sesh).oboe_rate_std];
		pop_rate_timbre.bass_rate = [nat_data(sesh).bass_rate];
		pop_rate_timbre.bass_rate_std = [nat_data(sesh).bass_rate_std];
		accuracy_all(iCF, irep) = accuracy;
	end
	%% Plot outputs
	% 
	% % Compute classification using all data
	% nexttile
	% C = confusionmat(pop_rate_timbre.response, pop_rate_timbre.predictions);
	% confusionchart(C)
	% title(sprintf('Accuracy = %0.2f%%, n=%d', pop_rate_timbre.accuracy*100,...
	% 	length(pop_rate_timbre.CFs)))
	% 
	% % Initial beta weight investigation
	% beta_weights = pop_rate_timbre.trainedClassifier.ClassificationSVM.Beta;
	% CFs = pop_rate_timbre.CFs;
	% 
	% nexttile
	% scatter(CFs/1000, abs(beta_weights), 'filled', 'MarkerEdgeColor','k', ...
	% 	'MarkerFaceAlpha',0.5)
	% hold on
	% yline(0)
	% xticks([0.1 0.2 0.5 1 2 5 10])
	% xlim([0.2 10])
	% xlabel('CF')
	% ylabel('abs(Beta Weight)')
	% set(gca, 'xscale', 'log')
	% title(CF_names{iCF})


end

%%

accuracies = [accuracy_all(2:7,:); accuracy_all(1,:)];
CF_names = {'Low', 'Medium', 'High', 'Low+Med', 'Low+High', 'Med+High','All'};

figure
boxplot(accuracies','Notch','on')
xlabel('Subset')
xticklabels(CF_names)
title('Model Accuracies (500 model repetitions)')
hold on

scattersize = 10;
for iCF = 1:7
	swarmchart(ones(size(accuracies, 2), 1)*iCF, accuracies(iCF,:), scattersize)
end

max_acc = max(accuracies, [], 2);
mean_acc = median(accuracies, 2);
plot(1:7, max_acc, 'k')
%plot(1:7, mean_acc, 'r')


% 
% [p,tbl,stats] = anova1(accuracy_all');
% results = multcompare(stats);

%% Save accuracies 

save(fullfile(base, 'model_comparisons', 'Pop_Rate_Timbre_Subset_CF.mat'), ...
	"accuracies", "CF_names")

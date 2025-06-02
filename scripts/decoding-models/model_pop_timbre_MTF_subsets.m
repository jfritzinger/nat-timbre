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
tiledlayout(5, 2)

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

% Get a subset of the data based on CF grouping
MTF_names = {'All', 'BE', 'BS', 'H', 'F'};

% Get a subset of the data based on MTF type
ind_all(1,:) = ones(1,num_data_all);
ind_all(2,:) = strcmp('BE', MTFs);
ind_all(3,:) = strcmp('BS', MTFs);
ind_all(4,:) = contains(MTFs,'H');
ind_all(5,:) = strcmp('F', MTFs);
ind_all = logical(ind_all);

for iMTF = 1:5

	for irep = 1:500
		ind = ind_all(iMTF,:);

		% Get random assortment of 60 units from each (for all, gets 20 per
		% group)
		if iMTF == 1
			ind_BE = find(ind_all(2,:));
			ind_BS = find(ind_all(3,:));
			ind_H = find(ind_all(4,:));
			ind_F = find(ind_all(5,:));

			rand_ind = [randsample(ind_BE, 6, false) ... % 51/246, 21% (6)
				randsample(ind_BS, 14, false)... % 123/246, 50% (14) 
				randsample(ind_H, 5, false)... % 43/246, 17% (5)
				randsample(ind_F, 3, false)]; % 29/246, 12% (3)
			%rand_ind = randsample(ind, 28,false);
		else
			rand_ind = randsample(find(ind), 28, false);
		end
		
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
		accuracy_all(iMTF, irep) = accuracy;
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
	% title(MTF_names{iMTF})


end

%%

figure
boxplot(accuracy_all','Notch','on')
xlabel('Subset')
xticklabels(MTF_names)
title('Model Accuracies (500 model repetitions)')
hold on
scattersize = 10;
for iCF = 1:5
	swarmchart(ones(size(accuracy_all, 2), 1)*iCF, accuracy_all(iCF,:), scattersize)
end
ylim([0.78 1])

% [p,tbl,stats] = anova1(accuracy_all');
% results = multcompare(stats);

%% Save accuracies 

save(fullfile(base, 'model_comparisons', 'Pop_Rate_Timbre_Subset_MTF.mat'), ...
	"accuracy_all", "MTF_names")

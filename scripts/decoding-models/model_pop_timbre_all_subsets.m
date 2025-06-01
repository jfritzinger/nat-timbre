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
tiledlayout(1, 4)

% Find all rows with bassoon and oboe
has_bass = ~cellfun(@isempty, {nat_data.bass_rate});
has_oboe = ~cellfun(@isempty, {nat_data.oboe_rate});
sesh_all = find(has_bass & has_oboe);
num_data_all = numel(sesh_all);
ind_b = 25:40;
ind_o = [1 3:17];
CFs_all1 = [nat_data(sesh_all).CF];
putative1 = {nat_data(sesh_all).putative};
MTFs1 = {nat_data(sesh_all).MTF};

ind_all(1,:) = strcmp('BE', MTFs1);
ind_all(2,:) = strcmp('BS', MTFs1);
ind_all(3,:) = contains(MTFs1,'H');
ind_all(4,:) = strcmp('F', MTFs1);
MTF_name = {'BE', 'BS', 'H', 'F'};

for iMTF = 1:4
	CFs_all = CFs_all1(ind_all(iMTF,:));
	putative = putative1(ind_all(iMTF,:));
	MTFs = MTFs1(ind_all(iMTF,:));

	% Get a subset of the data based on CF grouping
	CF_groups = [0, 14000; 0, 2000; 2000, 4000; 4000, 14000];
	CF_names = {'All', 'Low', 'Medium', 'High'};
	for iCF = 1:4

		for irep = 1:50
			ind = CFs_all > CF_groups(iCF, 1) & CFs_all < CF_groups(iCF, 2);

			% Subset!
			CFs = CFs_all(ind);
			sesh = sesh_all(ind);
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

	end

	%%

	nexttile
	boxplot(accuracy_all','Notch','on')
	xlabel('Subset')
	xticklabels(CF_names)
	title([MTF_name{iMTF} ' (50 reps, n=' num2str(length(ind)) ')'])
	hold on

	scattersize = 10;
	for iCF = 1:4
		swarmchart(ones(size(accuracy_all, 2), 1)*iCF, accuracy_all(iCF,:), scattersize)
	end
	ylim([0.64 1])
	%
	% [p,tbl,stats] = anova1(accuracy_all');
	% results = multcompare(stats);
end
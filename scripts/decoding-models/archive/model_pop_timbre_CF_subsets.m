%% model_pop_timbre_subsets
clear

%% Load in data

[base, ~, ~, ~] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')
load(fullfile(base, 'model_comparisons', 'Pop_Rate_Timbre_All.mat'), ...
	"pop_rate_timbre")

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

% Get best units
beta_weights = pop_rate_timbre.trainedClassifier.ClassificationSVM.Beta;
[~,originalpos] = sort(abs(beta_weights), 'descend' );
best_ind = originalpos(1:50);
[~,originalpos] = sort(abs(beta_weights), 'ascend' );
worst_ind = originalpos(1:50);

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

% % Get a subset of the data based on CF grouping
% MTF_names = {'All', 'BE', 'BS', 'H', 'F'};
% 
% % Get a subset of the data based on MTF type
% ind_all(1,:) = ones(1,num_data_all);
% ind_all(2,:) = strcmp('BE', MTFs);
% ind_all(3,:) = strcmp('BS', MTFs);
% ind_all(4,:) = contains(MTFs,'H');
% ind_all(5,:) = strcmp('F', MTFs);
% ind_all = logical(ind_all);

% Get a subset of the data based on CF grouping
% CF_groups = [0, 14000; 0, 1000; 1000 2000; 2000 3000;3000 4000; ...
% 	4000, 6000;6000 9000;9000 14000];
CF_groups = [0, 14000; 0, 2000; 2000 4000; 4000 8000];
%num_subset = 14;
%num_subset = 21;
%num_subset = 7;
num_subset = [1, 2, 3, 4, 5:5:90];
iii = 0;
for iCF = 1:length(CF_groups)
%for iMTF = 1:length(MTF_names)
	ind = CFs_all > CF_groups(iCF, 1) & CFs_all < CF_groups(iCF, 2);
	%ind = ind_all(iMTF,:);
	beta_weights = pop_rate_timbre.trainedClassifier.ClassificationSVM.Beta;
	[~,best_ind] = sort(abs(beta_weights(ind)), 'descend' );
	for inum = 1:length(num_subset)
		timerVal = tic;

		for irep = 1
			if num_subset(inum)>length(best_ind)
				accuracy_all(iCF, inum, irep) = NaN;
				%accuracy_all(iMTF, inum, irep) = NaN;
			else
				% Get random assortment of 60 units from each (for all, gets 20 per
				% group)
				% if iCF == 1
				% 	for itry = 1:7
				% 		ind_part{itry} = find(CFs_all > CF_groups(itry+1, 1) ...
				% 			& CFs_all < CF_groups(itry+1, 2));
				% 	end
				% 	rand_ind = [randsample(ind_part{1}, num_subset/7, false) ...
				% 		randsample(ind_part{2}, num_subset/7, false)...
				% 		randsample(ind_part{3}, num_subset/7, false)...
				% 		randsample(ind_part{4}, num_subset/7, false)...
				% 		randsample(ind_part{5}, num_subset/7, false)...
				% 		randsample(ind_part{6}, num_subset/7, false)...
				% 		randsample(ind_part{7}, num_subset/7, false)];
				% else % ismember(iCF, [2, 3, 4])
				% 	ind = CFs_all > CF_groups(iCF, 1) & CFs_all < CF_groups(iCF, 2);
				% 	%rand_ind = randsample(find(ind), num_subset, false);
				% % elseif iCF == 5
				% % 	rand_ind = [randsample(ind_low, num_subset/2, false) ...
				% % 		randsample(ind_med, num_subset/2, false)];
				% % elseif iCF == 6
				% % 	rand_ind = [randsample(ind_low, num_subset/2, false) ...
				% % 		randsample(ind_high, num_subset/2, false)];
				% % elseif iCF == 7
				% % 	rand_ind = [randsample(ind_med, num_subset/2, false) ...
				% % 		randsample(ind_high, num_subset/2, false)];
				% end
				
				%ind = find(CFs_all);
				%rand_ind = randsample(find(ind), num_subset(iCF), false);

				% Subset!
				rand_ind = best_ind(1:num_subset(inum));
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
				% pop_rate_timbre.trainedClassifier = trainedClassifier;
				% pop_rate_timbre.accuracy = accuracy;
				% pop_rate_timbre.predictions = predictions;
				% pop_rate_timbre.response = T.Instrument;
				% pop_rate_timbre.T = T;
				% pop_rate_timbre.CFs = CFs;
				% pop_rate_timbre.putative = putative;
				% pop_rate_timbre.sesh = sesh;
				% pop_rate_timbre.MTF = {nat_data(sesh).MTF};
				% pop_rate_timbre.oboe_rate = [nat_data(sesh).oboe_rate];
				% pop_rate_timbre.oboe_rate_std = [nat_data(sesh).oboe_rate_std];
				% pop_rate_timbre.bass_rate = [nat_data(sesh).bass_rate];
				% pop_rate_timbre.bass_rate_std = [nat_data(sesh).bass_rate_std];

				accuracy_all(iCF, inum, irep) = accuracy;
				%accuracy_all(iMTF, inum, irep) = accuracy;
			end
		end
		iii = iii+1;
		timer = toc(timerVal);
		fprintf('Models took %0.2g seconds, %d/%d, %0.0f%% done\n', timer, iii, 4*length(num_subset), iii/(4*length(num_subset))*100)
	end
end

%%

%accuracies = [accuracy_all(2:8,:); accuracy_all(1,:)];
% accuracies = accuracy_all;
% %CF_names = {'Low', 'Medium', 'High', 'Low+Med', 'Low+High', 'Med+High','All'};
% CF_names = {'<1000', '<2000', '<3000', '<4000', '<6000', '<8000', '<14000','All'};
% % CF_groups = [0, 14000; 0, 1000; 1000 2000; 2000 3000;3000 4000; ...
% % 	4000, 6000;6000 8000;8000 14000];
% figure
% boxplot(accuracies','Notch','on')
% xlabel('Subset')
% %xticklabels(CF_names)
% xticklabels(num_subset)
% title('Model Accuracies (500 model repetitions)')
% hold on
%
% scattersize = 10;
% for inum = 1:8
% 	swarmchart(ones(size(accuracies, 2), 1)*inum, accuracies(inum,:), scattersize)
% end
%
% max_acc = max(accuracies, [], 2);
% mean_acc = median(accuracies, 2);
% plot(1:8, max_acc, 'k')
% plot(1:8, mean_acc, 'r')

% [p,tbl,stats] = anova1(accuracy_all');
% results = multcompare(stats);

accuracy = mean(accuracy_all, 3);
std_acc = std(accuracy_all, [],3);

colorsMTF = {'#1b9e77', '#648FFF', '#DC267F', '#785EF0', '#FFB000'};
figure
hold on
% for iMTF = 1:5
% 	plot(num_subset, accuracy(iMTF,:), 'Color', colorsMTF{iMTF}, 'linewidth', 2)
% end
errorbar(num_subset, accuracy', std_acc'/sqrt(5), 'linewidth', 2)
ylabel('Accuracy')
xlabel('# Neurons in Model')
legend('All', 'CF = 0-2 kHz', 'CF = 2-4 kHz', 'CF = 4-8 kHz')
%legend(MTF_names)
grid on 
box off


%% Save accuracies

save(fullfile(base, 'model_comparisons', 'Pop_Rate_Timbre_Subset_CF2.mat'), ...
	"accuracy", "std_acc", "CF_groups")

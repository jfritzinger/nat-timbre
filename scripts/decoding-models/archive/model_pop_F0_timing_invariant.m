clear

%% Load in data
target = 'Bassoon';

base = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')
load(fullfile(base, 'model_comparisons', ['Neuron_Time_F0_' target '.mat']),...
	"neuron_time_F0")

%% Get correct output of model

tuning = readtable(fullfile(base, 'Tuning.xlsx')); % Load in tuning
target = 'Bassoon';
listing = dir(fullfile(base, 'waveforms', ['*' target '*.wav']));
files = {listing.name};
note_names = extractBetween(files, 'ff.', '.');
[~, index] = ismember(note_names, tuning.Note);
F0s1 = tuning.Frequency(index);
[F0s, order] = sort(F0s1);
% Get 180 matrix of both bassoon and timbre
has_bass = ~cellfun(@isempty, {nat_data.bass_rate});
has_oboe = ~cellfun(@isempty, {nat_data.oboe_rate});
sesh_all = find(has_bass & has_oboe);

% Find all rows with bassoon and oboe
ind_b = 25:40;
ind_o = [1 3:17];
F0s = F0s(ind_b);

putative_both = {nat_data(sesh_all).putative};
putative = {neuron_time_F0.putative};
for ii = 1:length(putative_both)
	put_target = putative_both{ii};
	index(ii) = find(cellfun(@(p) strcmp(put_target, p), putative));
end

% Get indices of interest
accuracy = [neuron_time_F0(index).accuracy];
[~,best_ind] = sort(abs(accuracy), 'descend' );
[~,worst_ind] = sort(abs(accuracy), 'ascend' );

% Find all rows with bassoon and oboe
num_data_all = numel(sesh_all);
CFs_all = [nat_data(sesh_all).CF];
putative = {nat_data(sesh_all).putative};
MTFs = {nat_data(sesh_all).MTF};


%% Get data

h_all2 = [];
%num_neurons = [1:4 5:5:40 50:10:100 150 200 1:4 5:5:40 50:10:100 150 200];
num_neurons = length(best_ind);
nmodels = length(num_neurons);
timerVal = tic;
for imodel = 1:nmodels
	timerVal = tic;

	for inrep = 1:1

		if imodel < nmodels/2+1 % 6 good models
			index = best_ind;
		else % 6 bad models
			index = worst_ind;
		end
		num_index = 1:num_neurons(imodel);
		num_data = num_neurons(imodel);
		sesh = sesh_all(index(num_index));

		T = getF0PopTable(nat_data, target, sesh, F0s, num_data, [], 'Timing');

		% Call model (classification)
		[trainedClassifier, validationAccuracy, validationPredictions] ...
			= trainClassifierPopTimeF0(T, F0s);

		% Plot
		%figure
		C = confusionmat(T.response, validationPredictions);
		%confusionchart(C)
		%title([target ', Accuracy = ' num2str(round(validationAccuracy*100)) '%'])

		% Calculate accuracy
		accuracy = sum(diag(C)) / sum(C(:)); % Calculate accuracy
		title(sprintf('%d neurons, Accuracy = %0.2f%%', num_data, accuracy*100))
		accur_all(imodel, inrep) = accuracy;
		C_all{imodel, inrep} = C;
	end
	timer = toc(timerVal);
	fprintf('Models took %0.2g minutes\n', timer/60)
	fprintf('%d/%d, %0.2f%% done!\n', imodel, nmodels, imodel/nmodels*100)
end

%%

% figure
% plot(num_neurons(1:nmodels/2), mean_acc(1:nmodels/2));
% hold on
% plot(num_neurons(nmodels/2+1:end), mean_acc(nmodels/2+1:end))
% xlabel('Number of Neurons in Model')
% ylabel('Accuracy')
% legend('Best', 'Worst')

mean_acc = mean(accur_all,2);
std_acc = std(accur_all, [], 2);

figure
errorbar(num_neurons(1:nmodels/2), mean_acc(1:nmodels/2), std_acc(1:nmodels/2)/sqrt(10));
hold on
errorbar(num_neurons(nmodels/2+1:end), mean_acc(nmodels/2+1:end), std_acc(nmodels/2+1:end)/sqrt(10))
xlabel('Number of Neurons in Model')
ylabel('Accuracy')
legend('Best', 'Worst')

%% Save data

save(fullfile(base, 'model_comparisons', 'pop_timing_F0_invariant_all.mat'), ...
	"accur_all","C_all", "num_neurons", "mean_acc", "std_acc", "trainedClassifier")
clear

%% Load in data
target = 'Bassoon';

base = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')
load(fullfile(base, 'model_comparisons', ['Neuron_Time_F0_' target '.mat']),...
	"neuron_time_F0")

%% Get correct output of model

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


% Get indices of interest
accuracy = [neuron_time_F0.accuracy];
[~,originalpos] = sort(abs(accuracy), 'descend' );
best_ind = originalpos(1:100);
[~,originalpos] = sort(abs(accuracy), 'ascend' );
worst_ind = originalpos(1:100);

%poolobj = parpool();
% figure
% bar(accuracy(originalpos))

%% Get data

h_all2 = [];
num_neurons = [1:4 5:5:50 1:4 5:5:50];
nmodels = length(num_neurons);
timerVal = tic;
for imodel = 1:nmodels
	timerVal = tic;


	if imodel < nmodels/2+1 % 6 good models
		index = best_ind;
	else % 6 bad models
		index = worst_ind;
	end
	num_index = 1:num_neurons(imodel);
	num_data = num_neurons(imodel);
	sesh = sesh_all(index(num_index));

	
	for ineuron = 1:num_data
		h_all = [];
		for itarget = 1:length(F0s)
			if strcmp(target, 'Bassoon')
				spikes_bass = nat_data(sesh(ineuron)).bass_spikerate{itarget}/1000; % ms
				spikereps_bass = nat_data(sesh(ineuron)).bass_spikerep{itarget};
			else
				spikes_bass = nat_data(sesh(ineuron)).oboe_spikerate{itarget}/1000; % ms
				spikereps_bass = nat_data(sesh(ineuron)).oboe_spikerep{itarget};
			end

			% Arrange data for SVM
			min_dis = 1;
			edges = 0:min_dis:300;
			t = 0+min_dis/2:min_dis:300-min_dis/2;
			for irep = 1:20
				h_bass(irep, :) = histcounts(spikes_bass(spikereps_bass==irep), edges);
			end
			h_all = [h_all; h_bass];
		end
		h_all2 = [h_all2, h_all];
	end

	% Put data into table
	response = reshape(repmat(F0s, 1, 20)', 1, []);
	T = array2table(h_all2);
	T.response = response';
	predictors = h_all2;

	% Call model (classification)
	[trainedClassifier, validationAccuracy, validationPredictions] ...
		= trainClassifierPopTimeF0(T, F0s);

	% Plot
	figure
	C = confusionmat(T.response, validationPredictions);
	confusionchart(C)
	title([target ', Accuracy = ' num2str(round(validationAccuracy*100)) '%'])

	% Calculate accuracy
	accuracy = sum(diag(C)) / sum(C(:)); % Calculate accuracy
	title(sprintf('%d neurons, Accuracy = %0.2f%%', num_data, accuracy*100))
	accur_all(imodel) = accuracy;

	timer = toc(timerVal);
	fprintf('Models took %0.2g minutes\n', timer/60)
	fprintf('%d/%d, %0.2f%% done!\n', imodel, nmodels, imodel/nmodels*100)	
end

%% 

figure
plot(num_neurons(1:nmodels/2), mean_acc(1:nmodels/2));
hold on
plot(num_neurons(nmodels/2+1:end), mean_acc(nmodels/2+1:end))
xlabel('Number of Neurons in Model')
ylabel('Accuracy')
legend('Best', 'Worst')

% mean_acc = mean(accuracies,2);
% std_acc = std(accuracies, [], 2);
%
% figure
% errorbar(num_neurons(1:nmodels/2), mean_acc(1:nmodels/2), std_acc(1:nmodels/2)/sqrt(10));
% hold on
% errorbar(num_neurons(nmodels/2+1:end), mean_acc(nmodels/2+1:end), std_acc(nmodels/2+1:end)/sqrt(10))
% xlabel('Number of Neurons in Model')
% ylabel('Accuracy')
% legend('Best', 'Worst')
clear

%% Load in data

base = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')
load(fullfile(base, 'model_comparisons', 'Neuron_Time_Timbre_All.mat'),...
	"neuron_time_timbre")

%% Get correct output of model

target = 'Bassoon';
tuning = readtable(fullfile(base, 'Tuning.xlsx')); % Load in tuning
listing = dir(fullfile(base, 'waveforms', '*.wav'));
target_WAV = arrayfun(@(n) contains(listing(n).name, target), 1:numel(listing), 'UniformOutput', false);
wav_nums =  find(cell2mat(target_WAV));

d = dir(fullfile(base,'waveforms', '*.wav'));
all_files = sort({d.name});
nfiles = length(wav_nums);
wav_npts = zeros(1,nfiles);
wav_data = cell(1,nfiles);
for i = 1:nfiles
	files{1,i} = all_files{wav_nums(i)};
end

% Sort by frequency of pitch
index = [];
note_names = extractBetween(files, 'ff.','.');
for ii = 1:nfiles % Find index of each note in tuning spreadsheet
	index(ii) = find(strcmp(note_names(ii), tuning.Note));
end
pitch_order = tuning.Frequency(index); % Get freqs of each note
[~, order] = sort(pitch_order); % Sort freqs
bass_pitch = pitch_order(order);


%% 

% Create array of correct responses
response = repmat([ones(1, 20) repmat(2, 1, 20)]', 16, 1);

% Find all rows with bassoon and oboe
has_bass = ~cellfun(@isempty, {nat_data.bass_rate});
has_oboe = ~cellfun(@isempty, {nat_data.oboe_rate});
sesh_all = find(has_bass & has_oboe);
num_data = numel(sesh_all);
ind_b = 25:40;
ind_o = [1 3:17];
CFs = [nat_data(sesh_all).CF];

% Get putative list
% putative = {nat_data(sesh_all).putative};
num_data_all = numel(sesh_all);
CFs_all = [nat_data(sesh_all).CF];
putative = {nat_data(sesh_all).putative};
MTFs = {nat_data(sesh_all).MTF};

% Get indices of interest
accuracy = [neuron_time_timbre.accuracy];
[~,originalpos] = sort(abs(accuracy), 'descend' );
best_ind = originalpos(1:100);
[~,originalpos] = sort(abs(accuracy), 'ascend' );
worst_ind = originalpos(1:100);


%% Get data

h_all2 = [];
num_neurons = [1:4 5:5:40 1:4 5:5:40];
nmodels = length(num_neurons);
timerVal = tic;
for imodel = 1:nmodels
	timerVal = tic;

	for inrep = 1:5

		if imodel < nmodels/2+1 % 6 good models
			index = best_ind;
		else % 6 bad models
			index = worst_ind;
		end
		num_index = 1:num_neurons(imodel);
		num_data = num_neurons(imodel);
		sesh = sesh_all(index(num_index));

		h_all2 = [];
		for ineuron = 1:num_data
			sesh_current = sesh(ineuron);
			h_all = [];
			h = [];
			for target = 1:16

				% Get data
				spikes_bass = nat_data(sesh_current).bass_spikerate{ind_b(target)}/1000; % ms
				spikereps_bass = nat_data(sesh_current).bass_spikerep{ind_b(target)};
				spikes_oboe = nat_data(sesh_current).oboe_spikerate{ind_o(target)}/1000;
				spikereps_oboe = nat_data(sesh_current).oboe_spikerep{ind_o(target)};

				% Arrange data for SVM
				min_dis = 1;
				edges = 0:min_dis:300;
				t = 0+min_dis/2:min_dis:300-min_dis/2;
				for irep = 1:20
					h_bass(irep, :) = histcounts(spikes_bass(spikereps_bass==irep), edges);
					h_oboe(irep, :) = histcounts(spikes_oboe(spikereps_oboe==irep), edges);
				end
				h_all = [h_bass; h_oboe];
				h = [h; h_all];
			end
			h_all2 = [h_all2, h];
		end

		% Put data into table
		T = array2table(h_all2);
		T.Instrument = response;
		predictors = h_all2;

		% Call model (classification)
		[trainedClassifier, validationAccuracy, validationPredictions] = ...
			trainClassifierNeuronTimeTimbre(T);

		% Plot
		%figure
		C = confusionmat(T.Instrument, validationPredictions);
		%confusionchart(C)

		% Calculate accuracy
		accuracy = sum(diag(C)) / sum(C(:)); % Calculate accuracy
		%title(sprintf('%d neurons, Accuracy = %0.2f%%', num_data, accuracy*100))
		accur_all(imodel, inrep) = accuracy;
		C_all{imodel, inrep} = C;
	end
	timer = toc(timerVal);
	fprintf('Models took %0.2g minutes\n', timer/60)
	fprintf('%d/%d, %0.2f%% done!\n', imodel, nmodels, imodel/nmodels*100)
end

%% Plot results 

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

save(fullfile(base, 'model_comparisons', 'pop_timing_timbre.mat'), ...
	"accur_all","C_all", "num_neurons", "mean_acc", "std_acc")


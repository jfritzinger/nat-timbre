%% model_neuron_rate_F0
clear


%% Load in data

[base, ~, ~, ~] = getPathsNT();
load(fullfile(base,'model_comparisons', 'Data_NT_3.mat'), 'nat_data')

%% Get list of all timbre stimuli (bassoon)

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


%% Shape data into model input

% Find all rows with bassoon in them
sesh = [];
for ii = 1:length(nat_data)
	rate = nat_data(ii).bass_VS;
	rate2 = nat_data(ii).oboe_VS;
	if ~isempty(rate) && ~isempty(rate2)
		sesh = [sesh ii];
	end
end
num_data = length(sesh);
ind_b = 25:40;
ind_o = [1 3:17];
target = 1;

%% Putting all data together

for ind = 1:num_data
	index = sesh(ind);

	h = [];
	for target = 1:16

		% Get data
		spikes_bass = nat_data(index).bass_spikerate{ind_b(target)}/1000; % ms
		spikereps_bass = nat_data(index).bass_spikerep{ind_b(target)};
		spikes_oboe = nat_data(index).oboe_spikerate{ind_o(target)}/1000;
		spikereps_oboe = nat_data(index).oboe_spikerep{ind_o(target)};

		% Arrange data for SVM
		min_dis = 300; %0.25;
		edges = 0:min_dis:300;
		t = 0+min_dis/2:min_dis:300-min_dis/2;
		for irep = 1:20
			h_bass1 = histcounts(spikes_bass(spikereps_bass==irep), edges);
			h_bass(irep, :) = h_bass1(randperm(length(h_bass1)));
			h_oboe1 = histcounts(spikes_oboe(spikereps_oboe==irep), edges);
			h_oboe(irep, :) = h_oboe1(randperm(length(h_bass1)));
		end
		h_all = [h_bass; h_oboe];
		h = [h; h_all];
	end

	% Put data into table
	T = array2table(h);
	T.Instrument = repmat([ones(20,1); ones(20, 1)*2], 16, 1);

	% Call SVM
	[trainedClassifier, validationAccuracy, validationPredictions] = ...
		trainClassifierTimeTimbre(T);
	C = confusionmat(validationPredictions, T.Instrument);
	accuracy_all(ind) = sum(diag(C)) / sum(C, "all");

	% Set up struct to save data
	neuron_time_timbre(ind).putative = nat_data(index).putative;
	neuron_time_timbre(ind).CF = nat_data(index).CF;
	neuron_time_timbre(ind).MTF = nat_data(index).MTF;
	neuron_time_timbre(ind).T = T;
	neuron_time_timbre(ind).Instrument = T.Instrument;
	neuron_time_timbre(ind).prediction = validationPredictions;
	neuron_time_timbre(ind).accuracy = accuracy_all(ind);
	neuron_time_timbre(ind).C = C;

	% h = [];
	% clear h_bass h_oboe
	% for target = 1:16
	% 
	% 	% Get data
	% 	spikes_bass = nat_data(index).bass_spikerate{ind_b(target)}/1000; % ms
	% 	spikereps_bass = nat_data(index).bass_spikerep{ind_b(target)};
	% 	spikes_oboe = nat_data(index).oboe_spikerate{ind_o(target)}/1000;
	% 	spikereps_oboe = nat_data(index).oboe_spikerep{ind_o(target)};
	% 
	% 	% Arrange data for SVM
	% 	min_dis = 5;
	% 	edges = 0:min_dis:300;
	% 	t = 0+min_dis/2:min_dis:300-min_dis/2;
	% 	for irep = 1:20
	% 		h_bass(irep, :) = histcounts(spikes_bass(spikereps_bass==irep), edges);
	% 		h_oboe(irep, :) = histcounts(spikes_oboe(spikereps_oboe==irep), edges);
	% 	end
	% 	h_all = [h_bass; h_oboe];
	% 	h = [h; h_all];
	% end
	% 
	% % Put data into table
	% T = array2table(h);
	% T.Instrument = repmat([ones(20,1); ones(20, 1)*2], 16, 1);
	% 
	% % Call SVM
	% [trainedClassifier, validationAccuracy, validationPredictions] = ...
	% 	trainClassifierTimeTimbre(T);
	% C = confusionmat(validationPredictions, T.Instrument);
	% accuracy_60(ind) = sum(diag(C)) / sum(C, "all");

	% Plots
	% figure
	% tiledlayout(2, 1)
	% nexttile
	% histogram(spikes_bass, 301)
	% nexttile
	% histogram(spikes_oboe, 301)
	% figure
	% tiledlayout(2, 1)
	% nexttile
	% hold on
	% for irep = 1:20
	% 	plot(t,h_bass(irep,:)+irep)
	% end
	% scatter(spikes_bass, spikereps_bass, 'filled')
	% nexttile
	% hold on
	% for irep = 1:20
	% 	plot(t,h_oboe(irep,:)+irep)
	% end
	% scatter(spikes_oboe, spikereps_oboe, 'filled')

	fprintf('%d/%d, %0.2f%% done! Accur=%0.2f%%\n', ind, ...
		num_data, ind/num_data*100, accuracy_all(ind)*100)
end

save(fullfile(base, 'model_comparisons', 'Neuron_Time_Timbre_All_Coarse300.mat'), ...
	"neuron_time_timbre", '-v7.3')

%% Plot accuracy of each neuron

figure
histogram(accuracy_all*100,21)
mean_F0 = mean(accuracy_all);
hold on
xline(50, 'k')
xline(mean_F0*100, 'r', 'LineWidth',2)
ylabel('# Neurons')
xlabel('Prediction Accuracy (%)')
%title(['F0=' num2str(round(bass_pitch(ind_b(ii))))])
xlim([0 100])
mean_all = mean(accuracy_all, 'all');
fprintf('Mean for all = %0.4f\n', mean_all)

%save('Neuron_Time_Timbre_All.mat', "accuracy", "mean_F0")

%% Plot scatter plot 

[base, datapath, savepath, ppi] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Neuron_Rate_Timbre_All.mat'), ...
	"neuron_rate_timbre")
accuracy_rate = [neuron_rate_timbre.accuracy];

figure
scatter(accuracy_rate, accuracy_all, 'filled', 'MarkerEdgeColor','k');
hold on
xlim([0.4 1])
ylim([0.4 1])
plot([0.4 1], [0.4 1], 'k')
ylabel('')
xlabel('Rate Accuracy')
ylabel('Timing Accuracy')
title('Model using shuffled              PSTH')

%% Get all rates for each repetition for bassoon (one example neuron)
%
% for ind = 1:num_data
% 	index = sesh(ind);
%
% 	for target = 1:16
%
% 		% Get data
% 		spikes_bass = nat_data(index).bass_spikerate{ind_b(target)}/1000; % ms
% 		spikereps_bass = nat_data(index).bass_spikerep{ind_b(target)};
% 		spikes_oboe = nat_data(index).oboe_spikerate{ind_o(target)}/1000;
% 		spikereps_oboe = nat_data(index).oboe_spikerep{ind_o(target)};
%
% 		% Arrange data for SVM
% 		min_dis = 0.25;
% 		edges = 0:min_dis:300;
% 		t = 0+min_dis/2:min_dis:300-min_dis/2;
% 		for irep = 1:20
% 			h_bass(irep, :) = histcounts(spikes_bass(spikereps_bass==irep), edges);
% 			h_oboe(irep, :) = histcounts(spikes_oboe(spikereps_oboe==irep), edges);
% 		end
% 		h_all = [h_bass; h_oboe];
%
% 		% Put data into table
% 		T = array2table(h_all);
% 		T.Instrument = [ones(20,1); ones(20, 1)*2];
%
% 		for imodelrep = 1
%
% 			% Call SVM
% 			SVMModel = fitcsvm(T,'Instrument','OptimizeHyperparameters','auto', ...
% 				'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName', ...
% 				'expected-improvement-plus', 'Verbose', 0));
% 			close gcf
% 			close gcf
%
% 			CVSVMModel = crossval(SVMModel, 'KFold', 5);
% 			classLoss = kfoldLoss(CVSVMModel);
%
% 			prediction = kfoldPredict(CVSVMModel);
% 			C = confusionmat(prediction, T.Instrument);
% 			accuracy(ind, target) = sum(diag(C)) / sum(C, "all");
% 		end
%
% 		% Set up struct to save data
% 		neuron_time_timbre(ind, target).putative = nat_data(index).putative;
% 		neuron_time_timbre(ind, target).ind_b = ind_b(target);
% 		neuron_time_timbre(ind, target).ind_o = ind_o(target);
% 		neuron_time_timbre(ind, target).F0 = bass_pitch(ind_b(target));
% 		neuron_time_timbre(ind, target).CF = nat_data(index).CF;
% 		neuron_time_timbre(ind, target).MTF = nat_data(index).MTF;
% 		neuron_time_timbre(ind, target).T = T;
% 		neuron_time_timbre(ind, target).Instrument = T.Instrument;
% 		neuron_time_timbre(ind, target).prediction = prediction;
% 		neuron_time_timbre(ind, target).accuracy = accuracy(ind,target);
% 		neuron_time_timbre(ind, target).C = C;
%
% 		% figure
% 		% confusionchart(C)
% 		% title(num2str(accuracy))
% 	end
% 	fprintf('%d/%d, %0.2f%% done!\n', ind, num_data, ind/num_data*100)
% end
%
% save(fullfile(base, 'model_comparisons', 'Neuron_Time_Timbre_Separate.mat'), ...
% 	"neuron_time_timbre", '-v7.3')
%
% % Plot accuracy of each neuron
% figure
% tiledlayout(4, 4)
% acc2 = max(accuracy, [], 3);
% for ii = 1:16
% 	nexttile
% 	histogram(acc2(:,ii)*100,21)
% 	mean_F0 = mean(acc2(:,ii));
% 	hold on
% 	xline(50, 'k')
% 	xline(mean_F0*100, 'r', 'LineWidth',2)
% 	ylabel('# Neurons')
% 	xlabel('Prediction Accuracy (%)')
% 	title(['F0=' num2str(round(bass_pitch(ind_b(ii))))])
% 	xlim([0 100])
% end
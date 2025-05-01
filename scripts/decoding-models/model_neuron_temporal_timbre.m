%% model_neuron_rate_F0
clear

%% Get list of all timbre stimuli (bassoon)

for iinstr = 2
	if iinstr == 1
		target = 'Oboe';
	else
		target = 'Bassoon';
	end

	if ismac
		fpath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data';
	else
		fpath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\';
	end
	tuning = readtable(fullfile(fpath, 'Tuning.xlsx')); % Load in tuning

	listing = dir(fullfile(fpath, 'waveforms', '*.wav'));
	target_WAV = arrayfun(@(n) contains(listing(n).name, target), 1:numel(listing), 'UniformOutput', false);
	wav_nums =  find(cell2mat(target_WAV));

	d = dir(fullfile(fpath,'waveforms', '*.wav'));
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

	if iinstr == 1
		files_o = files(order);
		note_names_o = note_names(order);
	else
		files_b = files(order);
		note_names_b = note_names(order);
	end
end
bass_pitch = pitch_order(order);

%% Load in data
if ismac
	filepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/model_comparisons';
else
	filepath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\model_comparisons';
end
load(fullfile(filepath, 'Data_NT.mat'), 'nat_data')

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

%% Get all rates for each repetition for bassoon (one example neuron)

for ind = 1:num_data
	index = sesh(ind);

	for target = 1:16

		% Get data 
		spikes_bass = nat_data(index).bass_spikerate{ind_b(target)}/1000; % ms
		spikereps_bass = nat_data(index).bass_spikerep{ind_b(target)};
		spikes_oboe = nat_data(index).oboe_spikerate{ind_o(target)}/1000;
		spikereps_oboe = nat_data(index).oboe_spikerep{ind_o(target)};

		% Arrange data for SVM 
		min_dis = 1;
		edges = 0:min_dis:300;
		t = 0+min_dis/2:min_dis:300-min_dis/2;
		for irep = 1:20
			h_bass(irep, :) = histcounts(spikes_bass(spikereps_bass==irep), edges);
			h_oboe(irep, :) = histcounts(spikes_oboe(spikereps_oboe==irep), edges);
		end
		h_all = [h_bass; h_oboe];

		% Put data into table 
		T = array2table(h_all);
		T.Instrument = [ones(20,1); ones(20, 1)*2];

		for imodelrep = 1:5
			% Call SVM
			SVMModel = fitcsvm(T,'Instrument','OptimizeHyperparameters','auto', ...
				'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName', ...
				'expected-improvement-plus', 'Verbose', 0));
			close gcf
			close gcf

			CVSVMModel = crossval(SVMModel, 'KFold', 5);
			classLoss = kfoldLoss(CVSVMModel);

			prediction = kfoldPredict(CVSVMModel);
			C = confusionmat(prediction, T.Instrument);
			accuracy(ind, target, imodelrep) = sum(diag(C)) / sum(C, "all");
		end
		% figure
		% confusionchart(C)
		% title(num2str(accuracy))

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


	end
	fprintf('%d/%d, %0.2f%% done!\n', ind, num_data, ind/num_data*100)
end

% Plot accuracy of each neuron
figure
tiledlayout(4, 4)
acc2 = max(accuracy, [], 3);
for ii = 1:16
	nexttile
	histogram(acc2(:,ii)*100,21)
	mean_F0 = mean(acc2(:,ii));
	hold on
	xline(50, 'k')
	xline(mean_F0*100, 'r', 'LineWidth',2)
	ylabel('# Neurons')
	xlabel('Prediction Accuracy (%)')
	title(['F0=' num2str(round(bass_pitch(ind_b(ii))))])
	xlim([0 100])
end

mean_all = mean(acc2, 'all');
fprintf('Mean for all = %0.4f\n', mean_all)

save('model_neuron_temporal_timbre.mat', "accuracy", "acc2", "mean_F0")

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

	% Put data into table
	T = array2table(h);
	T.Instrument = repmat([ones(20,1); ones(20, 1)*2], 16, 1);

	% Call SVM
	SVMModel = fitcsvm(T,'Instrument','OptimizeHyperparameters','auto', ...
		'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName', ...
		'expected-improvement-plus', 'Verbose', 0));

	CVSVMModel = crossval(SVMModel, 'KFold', 5);
	classLoss = kfoldLoss(CVSVMModel);

	prediction = kfoldPredict(CVSVMModel);
	C = confusionmat(prediction, T.Instrument);
	accuracy = sum(diag(C)) / sum(C, "all");
	accuracy_all(ind) = sum(diag(C)) / sum(C, "all");
	close all

	figure
	confusionchart(C)
	title(num2str(accuracy))

% 	% Plots
% 	% figure
% 	% tiledlayout(2, 1)
% 	% nexttile
% 	% histogram(spikes_bass, 301)
% 	% nexttile
% 	% histogram(spikes_oboe, 301)
% 	% figure
% 	% tiledlayout(2, 1)
% 	% nexttile
% 	% hold on
% 	% for irep = 1:20
% 	% 	plot(t,h_bass(irep,:)+irep)
% 	% end
% 	% scatter(spikes_bass, spikereps_bass, 'filled')
% 	% nexttile
% 	% hold on
% 	% for irep = 1:20
% 	% 	plot(t,h_oboe(irep,:)+irep)
% 	% end
% 	% scatter(spikes_oboe, spikereps_oboe, 'filled')

	fprintf('%d/%d, %0.2f%% done!\n', ind, num_data, ind/num_data*100)
end

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
save('model_neuron_temporal_timbre_all.mat', "accuracy", "mean_F0")


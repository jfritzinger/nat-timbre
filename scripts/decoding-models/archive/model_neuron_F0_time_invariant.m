%% model_neuron_time_F0_invariant
clear

%% Load in data

[base, datapath, savepath, ppi] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')


%% Get correct output of model
target = 'Invariant';

% Find all rows with bassoon and oboe
[sesh, num_data] = getF0Sessions(nat_data, 'Invariant');
CFs = [nat_data(sesh).CF];
putative = {nat_data(sesh).putative};

% Model including all F0s
F0s = getF0s('Invariant');
F0s = log10(F0s);
response = reshape(repmat(F0s, 1, 40)', 1, []);
response = response';
ind_b = 25:40;
ind_o = [1 3:17];

%% Get all timings for each repetition for bassoon (one example neuron)

for ind = 1:num_data
	index = sesh(ind);

	% Model including all F0s
	h_all = [];
	for itarget = 1:length(ind_b)

		spikes_bass = nat_data(index).bass_spikerate{ind_b(itarget)}/1000; % ms
		spikereps_bass = nat_data(index).bass_spikerep{ind_b(itarget)};
		spikes_oboe = nat_data(index).oboe_spikerate{ind_o(itarget)}/1000; % ms
		spikereps_oboe = nat_data(index).oboe_spikerep{ind_o(itarget)};

		% Arrange data for SVM
		min_dis = 1;
		edges = 0:min_dis:300;
		t = 0+min_dis/2:min_dis:300-min_dis/2;
		for irep = 1:20
			h_bass(irep, :) = histcounts(spikes_bass(spikereps_bass==irep), edges);
			h_oboe(irep, :) = histcounts(spikes_oboe(spikereps_oboe==irep), edges);
		end
		h_all = [h_all; h_bass; h_oboe];
	end
	T = array2table(h_all);
	T.response = response;

	[validationPredictions] = trainClassifierNeuronTimeF0(T, F0s);

	% Compute validation accuracy
	correctPredictions = (validationPredictions == response);
	isMissing = isnan(response);
	correctPredictions = correctPredictions(~isMissing);
	validationAccuracy(ind) = sum(correctPredictions)/length(correctPredictions);
	C = confusionmat(validationPredictions, response);

		% figure
		% confusionchart(C)
		% title(num2str(validationAccuracy*100))

	% Save data for each
	neuron_time_F0(ind).putative = nat_data(index).putative;
	neuron_time_F0(ind).CF = nat_data(index).CF;
	neuron_time_F0(ind).MTF = nat_data(index).MTF;
	neuron_time_F0(ind).response = response;
	neuron_time_F0(ind).T = T;
	neuron_time_F0(ind).validationPredictions = validationPredictions;
	neuron_time_F0(ind).accuracy = validationAccuracy(ind);
	neuron_time_F0(ind).C = C;

	fprintf('%d/%d, %0.2f%% done!\n', ind, num_data, ind/num_data*100)
end

%% Save struct of data

% save(fullfile(base, 'model_comparisons', ['Neuron_Time_F0_' target '.mat']), ...
% 	"neuron_time_F0", '-v7.3')

%% Plot 

accuracy_all = [neuron_time_F0(:).accuracy];
figure
nexttile
edges = linspace(0, 1, 51);
histogram(accuracy_all, edges)
ylabel('Number of Neurons')
xlabel('Accuracy')
title('Predicting F0 invariant of timbre')
hold on
xline(0.025)

[~,best_ind] = sort(abs(accuracy_all), 'ascend' );
for ii = 1:246
	C_acc(ii,:) = diag(neuron_time_F0(best_ind(ii)).C);
end
nexttile
pcolor(10.^F0s, 1:246, C_acc, 'EdgeColor','none', 'EdgeAlpha',0)
set(gca, 'xscale', 'log')
xlabel('F0 (Hz)')
ylabel('Neurons Sorted by Accuracy')
title('Accuracy')
colorbar

%% Load in data 

% Load in oboe
load(fullfile(base, 'model_comparisons', 'Neuron_Time_F0_Oboe.mat'), ...
	"neuron_time_F0")
neuron_time_F0_oboe = neuron_time_F0;
accuracy_time_oboe = [neuron_time_F0_oboe.accuracy];

% Load in bassoon
load(fullfile(base, 'model_comparisons', 'Neuron_Time_F0_Bassoon.mat'), ...
	"neuron_time_F0")
accuracy_time_bass = [neuron_time_F0.accuracy];

[~,best_ind] = sort(abs(accuracy_time_bass), 'ascend' );
for ii = 1:287
	C_acc = diag(neuron_time_F0(best_ind(ii)).C);
	C_all_sub(ii,:) = C_acc(ind_b);
end

figure
nexttile
pcolor(10.^F0s, 1:287, C_all_sub, 'EdgeColor','none', 'EdgeAlpha',0)
set(gca, 'xscale', 'log')
xticks([60 100 200 350 550])
xlabel('F0 (Hz)')
ylabel('Neurons Sorted by Accuracy')
title('Bassoon Accuracy',...
	'HorizontalAlignment','center')

%%

[~,best_ind] = sort(abs(accuracy_time_oboe), 'ascend' );
for ii = 1:length(best_ind)
	C_acc2 = diag(neuron_time_F0_oboe(best_ind(ii)).C);
	C_all_sub2(ii,:) = C_acc2(ind_o);
end
nexttile
pcolor(10.^F0s, 1:length(best_ind), C_all_sub2, 'EdgeColor','none', 'EdgeAlpha',0)
set(gca, 'xscale', 'log')
xticks([60 100 200 350 550])
xlabel('F0 (Hz)')
ylabel('Neurons Sorted by Accuracy')
title('Oboe Accuracy',...
	'HorizontalAlignment','center')

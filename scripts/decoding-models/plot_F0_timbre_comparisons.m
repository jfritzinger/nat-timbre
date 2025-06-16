%% plot_pitch_timbre_comparisons
clear

%%

[base, ~, ~, ~] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')
figure
tiledlayout(2, 2)
scattersize = 15;

%%  Single-unit rate

target = {'Bassoon', 'Oboe'};
for itarget = 1:2

	% Load data
	load(fullfile(base, 'model_comparisons', ['Neuron_Rate_F0_' ...
		target{itarget} '.mat']), "neuron_rate_F0")
	load(fullfile(base, 'model_comparisons', 'Neuron_Rate_Timbre_All.mat'), ...
		"neuron_rate_timbre")

	% Simplify neuron_rate_F0 to match timbre
	putatives_timbre = {neuron_rate_timbre.putative};
	putatives_F0 = {neuron_rate_F0.putative};

	for i = 1:numel(putatives_timbre)
		idx = find(strcmp(putatives_F0, putatives_timbre{i}));
		if ~isempty(idx)
			indices(i) = idx;
		else
			indices(i) = 0; % or NaN, to indicate not found
		end
	end
	indices(indices==0) = [];
	accuracy_timbre = [neuron_rate_timbre.accuracy];
	accuracy_F0 = [neuron_rate_F0(indices).accuracy];

	nexttile
	scatter(accuracy_timbre, accuracy_F0, scattersize, 'filled', 'MarkerEdgeColor','k');
	xlabel('Timbre Accuracy (rate)')
	ylabel([target{itarget} ' F0 Accuracy (rate)'])

	xlim([0 1])
	ylim([0 1])

end

%% Single-unit timing

target = {'Bassoon', 'Oboe'};
for itarget = 1:2

	% Load data
	load(fullfile(base, 'model_comparisons', ['Neuron_Time_F0_' ...
		target{itarget} '.mat']), "neuron_time_F0")
	load(fullfile(base, 'model_comparisons', 'Neuron_Time_Timbre_All.mat'), ...
		"neuron_time_timbre")

	% Simplify neuron_rate_F0 to match timbre
	putatives_timbre = {neuron_time_timbre.putative};
	putatives_F0 = {neuron_time_F0.putative};

	for i = 1:numel(putatives_timbre)
		idx = find(strcmp(putatives_F0, putatives_timbre{i}));
		if ~isempty(idx)
			indices(i) = idx;
		else
			indices(i) = 0; % or NaN, to indicate not found
		end
	end
	indices(indices==0) = [];
	accuracy_timbre = [neuron_time_timbre.accuracy];
	accuracy_F0 = [neuron_time_F0(indices).accuracy];

	nexttile
	scatter(accuracy_timbre, accuracy_F0, scattersize, 'filled', 'MarkerEdgeColor','k');
	xlabel('Timbre Accuracy (time)')
	ylabel([target{itarget} ' F0 Accuracy (time)'])
	xlim([0 1])
	ylim([0 1])
end

%% Population rate
% 
% target = {'Bassoon', 'Oboe'};
% for itarget = 1:2
% 
% 	% Load data
% 	load(fullfile(base, 'model_comparisons', ['Pop_Rate_F0_' ...
% 		target{itarget} '3.mat']), "pop_rate_F0")
% 	load(fullfile(base, 'model_comparisons', 'Pop_Rate_Timbre_All.mat'), ...
% 		"pop_rate_timbre")
% 
% 	% Simplify neuron_rate_F0 to match timbre
% 	putatives_timbre = pop_rate_timbre.putative;
% 	putatives_F0 = pop_rate_F0.putative;
% 
% 	for i = 1:numel(putatives_timbre)
% 		idx = find(strcmp(putatives_F0, putatives_timbre{i}));
% 		if ~isempty(idx)
% 			indices(i) = idx;
% 		else
% 			indices(i) = 0; % or NaN, to indicate not found
% 		end
% 	end
% 	indices(indices==0) = [];
% 	accuracy_timbre = [pop_rate_timbre.trainedClassifier.ClassificationSVM.Beta];
% 	accuracy_F0 = [pop_rate_F0.imp.ImportanceMean(indices)];
% 
% 	nexttile
% 	scatter(accuracy_timbre, accuracy_F0, 'filled', 'MarkerEdgeColor','k');
% 	xlabel('Pop Timbre Accuracy (time)')
% 	ylabel([target{itarget} ' Pop F0 Accuracy (time)'])
% 
% end
%% Model vs Data
clear
base = getPathsNT();

%% Load in Model_NT3

load(fullfile(base, 'model_comparisons',  'Model_NT3.mat'), 'nat_model')
load(fullfile(base, 'model_comparisons', 'Model_N_Rate_Timbre_All.mat'), ...
	"neuron_rate_timbre")
putative_model = {neuron_rate_timbre.putative};

%% Get overlapped indices and cut data down to model size
for ind = 1:length(putative_model)
	putative = putative_model{ind};
	index(ind) = find(cellfun(@(p) strcmp(p, putative), {nat_model.putative}));
end

for ind = 1:length(index)
	if isempty(nat_model(index(ind)).RM_R2)
		nat_model(index(ind)).RM_R2 = 0;
	end
end

RM_R2 = [nat_model(index).RM_R2]';
MTF_R2 = [nat_model(index).MTF_R2]';
both_R2 = mean([RM_R2, MTF_R2], 2);

good_MTF = MTF_R2 > 0.7;
%good_MTF = RM_R2 > 0.4;
%good_MTF = both_R2 > 0.7;

%% Neuron Timbre

load(fullfile(base, 'model_comparisons', 'Model_N_Rate_Timbre_All.mat'), ...
	"neuron_rate_timbre")
model_rate_timbre = neuron_rate_timbre;
accurate_rate_model = [model_rate_timbre.accuracy];
putative_model = {model_rate_timbre.putative};
MTFs = {model_rate_timbre.MTF};
isBE = strcmp(MTFs, 'BE');


load(fullfile(base, 'model_comparisons', 'Neuron_Rate_Timbre_All.mat'), ...
	"neuron_rate_timbre")
accurate_rate_data = [neuron_rate_timbre.accuracy];
putative_neuron = {neuron_rate_timbre.putative};

% Get overlapped indices and cut data down to model size
for ind = 1:length(putative_model)
	putative = putative_model{ind};
	index(ind) = find(cellfun(@(p) strcmp(p, putative), putative_neuron));
end
accurate_rate_data = accurate_rate_data(index);

% Plot
figure
tiledlayout(3, 2)
nexttile
hold on
% scatter(accurate_rate_data(~isBE), accurate_rate_model(~isBE), 'filled',...
% 	'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.5)
% scatter(accurate_rate_data(isBE), accurate_rate_model(isBE), 'filled',...
% 	'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.5)
scatter(accurate_rate_data(~good_MTF), accurate_rate_model(~good_MTF), 'filled',...
	'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.5, 'MarkerFaceColor',[0.4 0.4 0.4])
scatter(accurate_rate_data(good_MTF), accurate_rate_model(good_MTF), 'filled',...
	'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.8, 'MarkerFaceColor','r')


xlim([0.4 1])
ylim([0.4 1])
xlabel('Data Accuracy')
ylabel('Model Accuracy')
title('Timbre, Neuron Rate')
grid on
plot([0 1], [0 1], 'k')

mdl = fitlm(accurate_rate_data, accurate_rate_model);
p = mdl.Coefficients{2,4};
fprintf('Instrument, Rate: p=%0.4f\n', p)

[h,p,ci,stats] = ttest(accurate_rate_data, accurate_rate_model);
fprintf('Instrument, Rate: Dist, p=%0.4f\n', p)

%% Neuron Timbre: Timing

load(fullfile(base, 'model_comparisons', 'Model_N_Time_Timbre_All.mat'), ...
	"neuron_time_timbre")
model_time_timbre = neuron_time_timbre;
accurate_time_model = [model_time_timbre.accuracy];

load(fullfile(base, 'model_comparisons', 'Neuron_Time_Timbre_All.mat'), ...
	"neuron_time_timbre")
accurate_time_data = [neuron_time_timbre.accuracy];
accurate_time_data = accurate_time_data(index);

% Plot
nexttile
hold on
% scatter(accurate_time_data(~isBE), accurate_time_model(~isBE), 'filled',...
% 	'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.5)
% scatter(accurate_time_data(isBE), accurate_time_model(isBE), 'filled',...
% 	'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.5)

scatter(accurate_time_data(~good_MTF), accurate_time_model(~good_MTF), 'filled',...
	'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.5, 'MarkerFaceColor',[0.4 0.4 0.4])
scatter(accurate_time_data(good_MTF), accurate_time_model(good_MTF), 'filled',...
	'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.8, 'MarkerFaceColor','r')
xlim([0.4 1])
ylim([0.4 1])
xlabel('Data Accuracy')
ylabel('Model Accuracy')
title('Timbre, Neuron Timing')
grid on
plot([0 1], [0 1], 'k')
legend('BS', 'BE')

mdl = fitlm(accurate_time_data, accurate_time_model);
p = mdl.Coefficients{2,4};
fprintf('Instrument Identification, Timing: p=%0.4f\n', p)

%% Neuron F0: Rate

% figure
% tiledlayout(1, 2)
targets = {'Bassoon', 'Oboe'};
for itarget = 1:2
	target = targets{itarget};

	load(fullfile(base, 'model_comparisons', ['Model_N_Rate_F0_' target '.mat']), ...
		"neuron_rate_F0")
	model_rate_F0 = neuron_rate_F0;
	accurate_rate_F0_model = [model_rate_F0.accuracy];
	putative_model = {model_rate_F0.putative};
	MTFs = {model_rate_timbre.MTF};
	isBE = strcmp(MTFs, 'BE');

	load(fullfile(base, 'model_comparisons', ['Neuron_Rate_F0_' target '.mat']), ...
		"neuron_rate_F0")
	accuracy_rate_F0 = [neuron_rate_F0.accuracy];
	putative_neuron = {neuron_rate_F0.putative};

	% Get overlapped indices and cut data down to model size
	index = zeros(length(putative_model), 1);
	for ind = 1:length(putative_model)
		putative = putative_model{ind};
		index(ind) = find(cellfun(@(p) strcmp(p, putative), putative_neuron));
	end
	accuracy_rate_F0 = accuracy_rate_F0(index);

	nexttile
	hold on
	% scatter(accuracy_rate_F0(~isBE), accurate_rate_F0_model(~isBE), 'filled',...
	% 	'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.5)
	% scatter(accuracy_rate_F0(isBE), accurate_rate_F0_model(isBE), 'filled',...
	% 	'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.5)

	scatter(accuracy_rate_F0(~good_MTF), accurate_rate_F0_model(~good_MTF), 'filled',...
		'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.5, 'MarkerFaceColor',[0.4 0.4 0.4])
	scatter(accuracy_rate_F0(good_MTF), accurate_rate_F0_model(good_MTF), 'filled',...
		'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.8, 'MarkerFaceColor','r')
	xlim([0 0.5])
	ylim([0 0.5])
	xlabel('Data Accuracy')
	ylabel('Model Accuracy')
	title([target ' F0, Neuron Rate'])
	grid on
	plot([0 1], [0 1], 'k')

	mdl = fitlm(accuracy_rate_F0, accurate_rate_F0_model);
	p = mdl.Coefficients{2,4};
	fprintf('%s Identification, Rate: p=%0.4f\n', target, p)
end

%% Neuron F0: Timing

% figure
% tiledlayout(1, 2)
targets = {'Bassoon', 'Oboe'};
for itarget = 1:2
	target = targets{itarget};

	load(fullfile(base, 'model_comparisons', ['Model_Neuron_Time_F0_' target '.mat']), ...
		"neuron_time_F0")
	model_time_F0 = neuron_time_F0;
	accurate_time_F0_model = [model_time_F0.accuracy];
	putative_model = {model_time_F0.putative};
	MTFs = {model_rate_timbre.MTF};
	isBE = strcmp(MTFs, 'BE');

	load(fullfile(base, 'model_comparisons', ['Neuron_Time_F0_' target '.mat']), ...
		"neuron_time_F0")
	accuracy_time_F0 = [neuron_time_F0.accuracy];
	putative_neuron = {neuron_time_F0.putative};

	% Get overlapped indices and cut data down to model size
	index = zeros(length(putative_model), 1);
	for ind = 1:length(putative_model)
		putative = putative_model{ind};
		possible_ind = find(cellfun(@(p) strcmp(p, putative), putative_neuron));
		if ~isempty(possible_ind)
			index(ind) = possible_ind(1);
		else
			index(ind) = ind;
		end
	end
	accuracy_time_F0 = accuracy_time_F0(index);

	nexttile
	hold on
	% scatter(accuracy_time_F0(~isBE), accurate_time_F0_model(~isBE), ...
	% 	'filled', 'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.5)
	% scatter(accuracy_time_F0(isBE), accurate_time_F0_model(isBE), ...
	% 	'filled', 'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.5)

	scatter(accuracy_time_F0(~good_MTF), accurate_time_F0_model(~good_MTF), 'filled',...
		'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.5, 'MarkerFaceColor',[0.4 0.4 0.4])
	scatter(accuracy_time_F0(good_MTF), accurate_time_F0_model(good_MTF), 'filled',...
		'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.8, 'MarkerFaceColor','r')

	% Look for specific putative neurons
	array_both = [accuracy_time_F0; accurate_time_F0_model];
	mg_nb = find(array_both(1, :)<0.2 & array_both(2,:)>0.55);
	mg_ng = find(array_both(1, :)>0.5 & array_both(2,:)>0.5);

	mb_ng = find(array_both(1, :)>0.55 & array_both(2,:)<0.2);

	xlim([0 1])
	ylim([0 1])
	xlabel('Data Accuracy')
	ylabel('Model Accuracy')
	title([target ' F0, Neuron Timing'])
	grid on
	plot([0 1], [0 1], 'k')

	mdl = fitlm(accuracy_time_F0, accurate_time_F0_model);
	p = mdl.Coefficients{2,4};
	fprintf('%s Identification, Timing: p=%0.4f\n', target, p)

end


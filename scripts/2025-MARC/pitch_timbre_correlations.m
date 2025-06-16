% pitch_timbre_correlations
clear
save_fig = 1;

%%

[base, ~, ~, ppi] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')

%% Set up figure

figure('Position',[50 50 6.4*ppi, 6*ppi])
%tiledlayout(3, 4, 'TileIndexing','columnmajor')
linewidth = 2;
fontsize = 18;
titlesize = 20;
legsize = 20;
scattersize = 40;
spont_color = [0.4 0.4 0.4];
CF_color = [0.7 0.7 0.7];
tiledlayout(2, 2)

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

	h(itarget) = subplot(2, 2, itarget);
	scatter(accuracy_timbre, accuracy_F0, scattersize, 'filled', ...
		'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.5);
	%xlabel('Instrument Accuracy (Rate)')
	if itarget == 1
		ylabel('F0 Accuracy                                ')
	end
	hold on

	xlim([0 1])
	ylim([0 1])
	yticks(0:0.2:1)
	xticks(0:0.2:1)
	xticklabels([])
	if itarget == 2
		yticklabels([])
	end
	xtickangle(0)
	grid on

	mdl = fitlm(accuracy_timbre, accuracy_F0);
	x = linspace(0, 1, 20);
	y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
	plot(x, y, '--k', 'linewidth', 2)
	hleg = legend('', ...
		sprintf('p=%0.04f',mdl.Coefficients{2,4}), 'fontsize', legsize, ...
		'location', 'northwest', 'box', 'off');
	%hleg.ItemTokenSize = [12, 8];
	set(gca, 'fontsize', fontsize)

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

	h(itarget+2) = subplot(2, 2, itarget+2);
	scatter(accuracy_timbre, accuracy_F0, scattersize, 'filled', ...
		'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.5);
	if itarget == 1
		xlabel('                         Instrument Accuracy')
	end
	xlim([0 1])
	ylim([0 1])
	yticks(0:0.2:1)
	xticks(0:0.2:1)
	xtickangle(0)
	grid on
	hold on
	if itarget == 2
		yticklabels([])
	end

	mdl = fitlm(accuracy_timbre, accuracy_F0);
	x = linspace(0, 1, 20);
	y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
	plot(x, y, '--k', 'linewidth', 2)
	hleg = legend('', ...
		sprintf('p=%0.04f',mdl.Coefficients{2,4}), 'fontsize', legsize, ...
		'location', 'northwest', 'box', 'off');
	%hleg.ItemTokenSize = [12, 12];
	set(gca, 'fontsize', fontsize)
end

%% Arrange and annotate

left = [0.18 0.61];
bottom = linspace(0.13, 0.57, 2);
height = 0.35;
width = 0.35;

set(h(1), 'position', [left(1) bottom(2) width height])
set(h(2), 'position', [left(2) bottom(2) width height])
set(h(3), 'position', [left(1) bottom(1) width height])
set(h(4), 'position', [left(2) bottom(1) width height])

annotation('textbox',[0.26 0.94 0.25 0.058],...
	'String','Bassoon','FontSize',titlesize,...
	'EdgeColor','none', 'FontWeight','bold');
annotation('textbox',[0.72 0.94 0.25 0.058],...
	'String','Oboe','FontSize',titlesize,...
	'EdgeColor','none', 'FontWeight','bold');
annotation('textbox',[0.05 0.7 0.25 0.058],...
	'String','Rate','FontSize',titlesize,...
	'EdgeColor','none', 'Rotation',90, 'FontWeight','bold');
annotation('textbox',[0.05 0.2 0.25 0.058],...
	'String','Timing','FontSize',titlesize,...
	'EdgeColor','none', 'Rotation',90, 'FontWeight','bold');


%% Save figure

if save_fig == 1
	filename = 'pitch_timbre_corelations';
	save_figure_MARC(filename)
end

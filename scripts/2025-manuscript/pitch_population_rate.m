%% plot_pitch_population_rate
clear
save_fig = 0;

%% Set up figure
[base, ~, ~, ppi] = getPathsNT();

figure('Position',[50,50,6.8*ppi,4.5*ppi])
%tiledlayout(2, 6, 'TileSpacing','tight', 'Padding','compact')
linewidth = 1;
fontsize = 8;
legsize = 7;
titlesize = 10;
labelsize = 12;
scattersize = 10;
colorsMTF = {'#648FFF', '#DC267F', '#785EF0', '#FFB000'};

%% Plot each row as bassoon or oboe

targets = {'Bassoon', 'Oboe', 'Invariant'};
ind = [0, 5, 10];

for iinstru = 1:3
	target = targets{iinstru};

	%% Load in data

	load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')
	load(fullfile(base, 'model_comparisons', ['Pop_Rate_F0_' target '.mat']),...
		"pop_rate_F0")
	load(fullfile(base, 'model_comparisons', ['Pop_Rate_F0_Linear_' target '.mat']),...
		"pred_F0", "T", "r2", "C", "accuracy", "Mdl")
	load(fullfile(base, 'model_comparisons', ['Pop_Rate_F0_Best_' target '.mat']), ...
		"accuracies", "num_neurons")

	%% A. Plot confusion matrix

	% Get stimulus
	tuning = readtable(fullfile(base, 'Tuning.xlsx')); % Load in tuning
	listing = dir(fullfile(base, 'waveforms', ['*' target '*.wav']));
	files = {listing.name};
	note_names = extractBetween(files, 'ff.', '.');
	[~, index] = ismember(note_names, tuning.Note);
	F0s1 = tuning.Frequency(index);
	[F0s, order] = sort(F0s1);
	if strcmp(target, 'Invariant')
		listing = dir(fullfile(base, 'waveforms', ['*' 'Bassoon' '*.wav']));
		files = {listing.name};
		note_names = extractBetween(files, 'ff.', '.');
		[~, index] = ismember(note_names, tuning.Note);
		F0s1 = tuning.Frequency(index);
		[F0s, order] = sort(F0s1);
		F0s = F0s(25:40);
	end

	h(ind(iinstru)+1) = subplot(3, 5, ind(iinstru)+1);
	pcolor(F0s, F0s, pop_rate_F0.C, 'EdgeColor','none')
	set(gca, 'xscale', 'log', 'yscale', 'log')
	clim([0 20])
	ylabel('Actual F0 (Hz)')
	if iinstru == 3
		xlabel('Predicted F0 (Hz)')
	end
	a=colorbar;
	a.Label.String = '# Accurate Predictions';
	xticks([50 100 200 500 1000 1600])
	xticklabels([50 100 200 500 1000 1600]/1000)
	yticks([50 100 200 500 1000 1600])
	yticklabels([50 100 200 500 1000 1600]/1000)
	set(gca, 'FontSize', fontsize)
	clim([0 10])
	colormap(brewermap([],"Blues"))

	%% B. Plot weights

	beta_mean = pop_rate_F0.imp.ImportanceMean;
	beta_std = pop_rate_F0.imp.ImportanceStandardDeviation;

	[~,originalpos] = sort(abs(beta_mean), 'descend' );
	h(ind(iinstru)+2) = subplot(3, 5, ind(iinstru)+2);
	ordered_mean = beta_mean(originalpos);
	ordered_std = beta_std(originalpos);

	bar(ordered_mean)
	hold on
	errorbar(1:length(beta_mean), ordered_mean, ordered_std/sqrt(10), ...
		'CapSize',2, 'color', 'k', 'LineStyle', 'none');
	xlim([0 25])
	box off
	set(gca, 'FontSize', fontsize)
	if iinstru == 3
		xlabel('Neuron #')
	end
	ylabel('Importance')
	grid on

	%% C. Plot CFs

	h(ind(iinstru)+3) = subplot(3, 5, ind(iinstru)+3);
	beta_mean = pop_rate_F0.imp.ImportanceMean;
	CFs = [pop_rate_F0.CF];
	MTFs = pop_rate_F0.MTF;
	isMTF(1,:) = strcmp(MTFs, 'BE');
	isMTF(2,:) = strcmp(MTFs, 'BS');
	isMTF(3,:) = contains(MTFs, 'H');
	isMTF(4,:) = strcmp(MTFs, 'F');

	hold on
	for iMTF = 1:4
		scatter(CFs(isMTF(iMTF, :)), beta_mean(isMTF(iMTF, :)), 10, 'filled', ...
			'MarkerEdgeColor','k', 'MarkerFaceColor',colorsMTF{iMTF})
	end
	set(gca, 'xscale', 'log')
	xticks([100 200 500 1000 2000 5000 10000])
	xticklabels([100 200 500 1000 2000 5000 10000]/1000)
	yticklabels([])
	set(gca, 'fontsize', fontsize)
	grid on
	if iinstru == 1
		legend('BE', 'BS', 'H', 'F', 'fontsize', legsize)
	elseif iinstru == 3
		xlabel('CF (kHz)')
	end

	clear isMTF
	set(gca, 'FontSize', fontsize)

	%% D. Plot subset
	mean_acc = mean(accuracies,2);
	std_acc = std(accuracies, [], 2);
	nmodels = length(num_neurons);

	h(ind(iinstru)+4) = subplot(3, 5, ind(iinstru)+4);
	errorbar(num_neurons(1:nmodels/2), mean_acc(1:nmodels/2), std_acc(1:nmodels/2)/sqrt(10));
	hold on
	errorbar(num_neurons(nmodels/2+1:end), mean_acc(nmodels/2+1:end), std_acc(nmodels/2+1:end)/sqrt(10))
	ylabel('Accuracy')
	if iinstru == 1
		hleg = legend('Best Neurons', 'Worst Neurons', 'fontsize', legsize);
		hleg.ItemTokenSize = [8, 8];
	elseif iinstru == 3
		xlabel('# Neurons')
	end
	set(gca, 'FontSize', fontsize)
	box off
	grid on

	%% E. Plot linear

	h(ind(iinstru)+5) = subplot(3, 5, ind(iinstru)+5);
	scatter(10.^(T.response), 10.^(pred_F0), scattersize, 'filled', ...
		'MarkerFaceAlpha',0.3, 'MarkerFaceColor','k')
	set(gca, 'xscale', 'log', 'yscale', 'log')
	hold on
	plot(10.^T.response, 10.^T.response, 'k')

	mdl = fitlm(T.response, pred_F0);
	x = linspace(58, 1661, 20);
	p(1) = mdl.Coefficients.Estimate(2,1);
	p(2) = mdl.Coefficients.Estimate(1,1);
	p(3) = mdl.Coefficients.pValue(2);
	y = 10^p(2) .* x.^p(1);
	%y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
	plot(x, y, 'r')

	if iinstru == 3
		xlabel('Actual F0 (Hz)')
	end
	ylabel('Predicted F0 (Hz)')
	set(gca, 'FontSize', fontsize)
	% hleg = legend('Fits', 'Unity', 'Fit Line');
	% hleg.ItemTokenSize = [8, 8];
	xlim([F0s(1) F0s(end)])
	xticks([10 20 50 100 200 500 1000 1600])
	xticklabels([10 20 50 100 200 500 1000 1600]/1000)
	yticks([10 20 50 100 200 500 1000 1600])
	yticklabels([10 20 50 100 200 500 1000 1600]/1000)
	grid on

	% Annotate
	msg = ['R^2 = ' num2str(round(r2, 2))];
	text(0.05, 0.95, msg, 'Units', 'normalized', ...
		'VerticalAlignment', 'top', 'FontSize',legsize)

	% %% F.
	%
	% closest_cat = [];
	% for ii = 1:length(pred_F0)
	% 	differences = abs(F0s - pred_F0(ii));
	% 	[~, closest_column_index] = min(differences);
	% 	closest_cat(ii) = F0s(closest_column_index);
	% end
	%
	% % Plot confusion matrix
	% h(ind(iinstru)+6) = subplot(2, 6, ind(iinstru)+6);
	% imagesc(C)
	% %title(sprintf('Accuracy = %0.2f%%', accuracy*100))
	% set(gca, 'FontSize', fontsize)
	% if iinstru == 2
	% 	xlabel('Predicted F0 (Hz)')
	% end
	% ylabel('Actual F0 (Hz)')
	% xticks([50 100 200 500 1000 1600])
	% xticklabels([50 100 200 500 1000 1600]/1000)
	% yticks([50 100 200 500 1000 1600])
	% yticklabels([50 100 200 500 1000 1600]/1000)
end

%% Arrange plots

left = [0.09 0.37 0.48 0.64 0.86];
bottom = fliplr(linspace(0.1, 0.74, 3));
width = 0.12;
height = 0.22;

for ii = 1:3
	set(h(ind(ii)+1), 'position', [left(1) bottom(ii) 0.14 height]);
	set(h(ind(ii)+2), 'position', [left(2) bottom(ii) 0.1 height]);
	set(h(ind(ii)+3), 'position', [left(3) bottom(ii) 0.1 height]);
	set(h(ind(ii)+4), 'position', [left(4) bottom(ii) width height]);
	set(h(ind(ii)+5), 'position', [left(5) bottom(ii) width height]);
	%set(h(ind(ii)+6), 'position', [left(6) bottom(ii) width height]);
end

%% Annotate

annotation("textbox", [0.07 0.75 0.1386 0.1088], "String", "Bassoon",...
	"FontSize", titlesize+2, "FontWeight", "bold", "EdgeColor", "none", "Rotation",90)
annotation("textbox", [0.03 0.46 0.09778 0.05666], "String", "Oboe",...
	"FontSize", titlesize+2, "FontWeight", "bold", "EdgeColor", "none", "Rotation",90)
annotation("textbox", [0.07 0.11 0.1386 0.1088], "String", "Invariant",...
	"FontSize", titlesize+2, "FontWeight", "bold", "EdgeColor", "none", "Rotation",90)

labelleft= left-0.055;
labelbottom = 0.94;
annotation('textbox',[labelleft(1) labelbottom 0.071 0.058],...
	'String','A','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(2) labelbottom 0.071 0.058],...
	'String','B','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(4) labelbottom 0.071 0.058],...
	'String','C','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(5) labelbottom 0.071 0.058],...
	'String','D','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');

%% Save figure

if save_fig == 1
	filename = 'fig7_pitch_pop_rate';
	save_figure(filename)
end

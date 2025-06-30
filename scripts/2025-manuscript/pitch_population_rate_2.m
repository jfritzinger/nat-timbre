%% plot_pitch_population_rate
clear
save_fig = 0;

%% Set up figure
[base, ~, ~, ppi] = getPathsNT();

figure('Position',[50,50,6.8*ppi,4*ppi])
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
	load(fullfile(base, 'model_comparisons', ['Pop_Rate_F0_' target '_CF.mat']), ...
		"accuracy", "num_subset")
	accuracies_CF = accuracy;
	load(fullfile(base, 'model_comparisons', ['Pop_Rate_F0_' target '_MTF.mat']), ...
		"accuracy", "num_subset", "MTF_names")
	accuracies = accuracy;


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

	%% B. Plot linear

	h(ind(iinstru)+2) = subplot(3, 5, ind(iinstru)+2);
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

	%% C. Plot CFs
	F0s = getF0s(target);
	F0s = log10(F0s);
	[sesh, num_data] = getF0Sessions(nat_data, target);
	T = getF0PopTable(nat_data, target, sesh, F0s, num_data, 'linear');
	CFs = [nat_data(sesh).CF];
	MTFs = {nat_data(sesh).MTF};
	putative = {nat_data(sesh).putative};

	trainedModels = Mdl.Trained;  % Cell array of models (one per fold)

	% Preallocate matrix for coefficients
	numFeatures = size(trainedModels{1}.Beta, 1);
	allBetas = zeros(numFeatures, Mdl.KFold);

	% Populate coefficients from each fold
	for fold = 1:Mdl.KFold
		allBetas(:, fold) = trainedModels{fold}.Beta;
	end
	avgBeta = mean(allBetas, 2);  % Average across folds
	h(ind(iinstru)+3) = subplot(3, 5, ind(iinstru)+3);

	scatter(CFs, avgBeta, scattersize, 'filled', 'MarkerEdgeColor','k');
	hold on
	set(gca, 'xscale', 'log')

	mdl = fitlm(CFs, avgBeta);
	x = linspace(300, 14000, 50);
	y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
	plot(x, y, ':k')
	hleg = legend('', ...
		sprintf('p=%0.04f',mdl.Coefficients{2,4}), 'fontsize', legsize, ...
		'location', 'northwest', 'box', 'off');
	hleg.ItemTokenSize = [8, 8];
	if iinstru ==3
		xlabel('CFs')
	end
	ylabel('Beta Weights')
	xticks([200 500 1000 2000 5000 10000])
	xticklabels([0.2 0.5 1 2 5 10])
	grid on
	set(gca, 'FontSize', fontsize)

	%% D. Plot subset CF
	mean_acc = accuracies_CF;
	nmodels = size(accuracies_CF, 2);

	h(ind(iinstru)+4) = subplot(3, 5, ind(iinstru)+4);
	plot(num_subset(1:nmodels/2), mean_acc(:, 1:nmodels/2));
	hold on
	%errorbar(num_subset(nmodels/2+1:end), mean_acc(nmodels/2+1:end), std_acc(nmodels/2+1:end)/sqrt(10))
	ylabel('R^2')
	if iinstru == 1
		hleg = legend('Best', 'Low CF', 'Med CF', 'High CF', ...
			'fontsize', legsize, 'box', 'off');
		hleg.ItemTokenSize = [8, 8];
	elseif iinstru == 3
		xlabel('# Neurons')
	end
	set(gca, 'FontSize', fontsize)
	box off
	grid on

	%% E. Plot subset MTF
	mean_acc = accuracies;
	nmodels = size(accuracies, 2);

	h(ind(iinstru)+5) = subplot(3, 5, ind(iinstru)+5);
	plot(num_subset(1:nmodels/2), mean_acc(:, 1:nmodels/2));
	hold on
	%errorbar(num_neurons(nmodels/2+1:end), mean_acc(nmodels/2+1:end), std_acc(nmodels/2+1:end)/sqrt(10))
	if iinstru == 1
		hleg = legend(['Best', MTF_names(2:end)], 'fontsize', ...
			legsize, 'box', 'off');
		hleg.ItemTokenSize = [8, 8];
	elseif iinstru == 3
		xlabel('# Neurons')
	end
	set(gca, 'FontSize', fontsize)
	box off
	grid on

end

%% Arrange plots

left = [ 0.0900    0.3700    0.55    0.72    0.8700];
bottom = fliplr(linspace(0.09, 0.71, 3));
width = 0.12;
height = 0.22;

for ii = 1:3
	set(h(ind(ii)+1), 'position', [left(1) bottom(ii) 0.14 height]);
	set(h(ind(ii)+2), 'position', [left(2) bottom(ii) 0.1 height]);
	set(h(ind(ii)+3), 'position', [left(3) bottom(ii) 0.1 height]);
	set(h(ind(ii)+4), 'position', [left(4) bottom(ii) width height]);
	set(h(ind(ii)+5), 'position', [left(5) bottom(ii) width height]);
end

%% Annotate

annotation("textbox", [0.06 0.72 0.1386 0.1088], "String", "Bassoon",...
	"FontSize", titlesize+2, "FontWeight", "bold", "EdgeColor", "none", "Rotation",90)
annotation("textbox", [0.03 0.44 0.09778 0.05666], "String", "Oboe",...
	"FontSize", titlesize+2, "FontWeight", "bold", "EdgeColor", "none", "Rotation",90)
annotation("textbox", [0.06 0.09 0.1386 0.1088], "String", "Invariant",...
	"FontSize", titlesize+2, "FontWeight", "bold", "EdgeColor", "none", "Rotation",90)

labelleft= left-0.05;
labelbottom = 0.95;
annotation('textbox',[labelleft(1) labelbottom 0.071 0.058],...
	'String','A','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(2) labelbottom 0.071 0.058],...
	'String','B','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(3) labelbottom 0.071 0.058],...
	'String','C','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(4) labelbottom 0.071 0.058],...
	'String','D','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(5) labelbottom 0.071 0.058],...
	'String','E','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');

%% Save figure

if save_fig == 1
	filename = 'fig7_pitch_pop_rate';
	save_figure(filename)
end

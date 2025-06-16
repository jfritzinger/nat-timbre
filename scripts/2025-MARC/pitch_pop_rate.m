%% pitch_pop_timing 
clear
save_fig = 1;

%% Set up figure

[base, ~, ~, ppi] = getPathsNT();
figure('Position',[50, 50, 12.5*ppi, 5*ppi])
linewidth = 2;
fontsize = 18;
titlesize = 20;
legsize = 20;
scattersize = 40;

%% Plot each row as bassoon or oboe

targets = {'Bassoon', 'Oboe'};
ind = [0, 3];

for iinstru = 1:2
	target = targets{iinstru};

	%% Load in data

	load(fullfile(base, 'model_comparisons', ['Pop_Rate_F0_' target '3.mat']),...
		"pop_rate_F0")
	load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')
	load(fullfile(base, 'model_comparisons', ['Pop_Rate_F0_Best_' target '.mat']), ...
		"accuracies", "num_neurons")
	load(fullfile(base, 'model_comparisons', ['pop_rate_F0_linear_' target '.mat']),...
		"pred_F0", "T_test", "T", "r2", "C", "accuracy", "F0s")

	%% A. Plot confusion matrix

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
	
	h(ind(iinstru)+1) = subplot(2, 3, ind(iinstru)+1);
	pcolor(F0s, F0s, pop_rate_F0.C, 'EdgeColor','none')
	set(gca, 'xscale', 'log', 'yscale', 'log')
	clim([0 10])
	if iinstru == 1
		ylabel('Actual F0 (Hz)                           ')
	end
	if iinstru == 2
		xlabel('Predicted F0 (Hz)')
	end
	colormap(brewermap([],"Blues"))

	a=colorbar;
	if iinstru == 1
		a.Label.String = '# Accurate Predictions                                           ';
	end
	xticks([50 100 200 500 1000 1600])
	xticklabels([50 100 200 500 1000 1600]/1000)
	yticks([50 100 200 500 1000 1600])
	yticklabels([50 100 200 500 1000 1600]/1000)
	set(gca, 'FontSize', fontsize)


	%% D. Plot subset
	mean_acc = mean(accuracies,2);
	std_acc = std(accuracies, [], 2);
	nmodels = length(num_neurons);

	h(ind(iinstru)+2) = subplot(2, 3, ind(iinstru)+2);
	errorbar(num_neurons(1:nmodels/2), mean_acc(1:nmodels/2), ...
		std_acc(1:nmodels/2)/sqrt(10), 'LineWidth',2);
	hold on
	if iinstru == 1
		ylabel('Accuracy                              ')
	end
	if iinstru == 1
		hleg = legend('Best Neurons', 'fontsize', legsize, 'location', 'southeast');
		hleg.ItemTokenSize = [8, 8];
	else
		xlabel('# Neurons')
	end
	set(gca, 'FontSize', fontsize)
	box off
	grid on

	%% E. Plot linear

	h(ind(iinstru)+3) = subplot(2, 3, ind(iinstru)+3);
	scatter(10.^(T.response), 10.^(pred_F0), scattersize, 'filled', 'MarkerFaceAlpha',0.3, 'MarkerFaceColor','k')
	set(gca, 'xscale', 'log', 'yscale', 'log')
	hold on
	plot(10.^T_test.response, 10.^T_test.response, 'k')

	mdl = fitlm(T.response, pred_F0);
	x = linspace(58, 1661, 20);
	p(1) = mdl.Coefficients.Estimate(2,1);
	p(2) = mdl.Coefficients.Estimate(1,1);
	p(3) = mdl.Coefficients.pValue(2);
	y = 10^p(2) .* x.^p(1);
	%y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
	% plot(x, y, 'r')

	if iinstru == 2
		xlabel('Actual F0 (kHz)')
	end
	if iinstru == 1
		ylabel('Predicted F0 (kHz)                        ')
	end
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

end

%% Arrange plots

left = [0.12 0.47 0.78];
bottom = [0.59 0.14];
width = 0.2;
height = 0.36;

for ii = 1:2
	set(h(ind(ii)+1), 'position', [left(1) bottom(ii) width height]);
	set(h(ind(ii)+2), 'position', [left(2) bottom(ii) width height]);
	set(h(ind(ii)+3), 'position', [left(3) bottom(ii) width height]);
end

%% Annotate

annotation("textbox", [0.05 0.62 0.1386 0.1088], "String", "Bassoon",...
	"FontSize", titlesize+2, "FontWeight", "bold", "EdgeColor", "none", "Rotation",90)
annotation("textbox", [0.025 0.22 0.09778 0.05666], "String", "Oboe",...
	"FontSize", titlesize+2, "FontWeight", "bold", "EdgeColor", "none", "Rotation",90)

%% Save figure

if save_fig == 1
	filename = 'pitch_pop_rate';
	save_figure_MARC(filename)
end

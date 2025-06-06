%% plot_pitch_population_rate
clear

%% Set up figure
[base, ~, ~, ppi] = getPathsNT();

figure('Position',[50,50,10*ppi,3.5*ppi])
tiledlayout(2, 6, 'TileSpacing','tight', 'Padding','compact')
linewidth = 1;
fontsize = 8;
legsize = 7;
labelsize = 12;
scattersize = 10;

%% Plot each row as bassoon or oboe

targets = {'Bassoon', 'Oboe'};

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
	nexttile
	imagesc(pop_rate_F0.C)
	xlabel('Predicted Class')
	ylabel('Actual Class')
	set(gca, 'FontSize', fontsize)
	title(['Accuracy = ' num2str(pop_rate_F0.validationAccuracy*100) '%'])

	%% B. Plot weights

	beta_mean = pop_rate_F0.imp.ImportanceMean;
	beta_std = pop_rate_F0.imp.ImportanceStandardDeviation;

	[~,originalpos] = sort(abs(beta_mean), 'descend' );
	nexttile
	ordered_mean = beta_mean(originalpos);
	ordered_std = beta_std(originalpos);

	bar(ordered_mean)
	hold on
	errorbar(1:length(beta_mean), ordered_mean, ordered_std/sqrt(10), ...
		'CapSize',2, 'color', 'k', 'LineStyle', 'none');
	xlim([0 100])
	box off
	set(gca, 'FontSize', fontsize)
	xlabel('Neuron #')
	ylabel('Importance')
	grid on

	%% C. Plot CFs

	nexttile
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
			'MarkerEdgeColor','k')
	end
	set(gca, 'xscale', 'log')
	xticks([100 200 500 1000 2000 5000 10000])
	xticklabels([100 200 500 1000 2000 5000 10000]/1000)
	ylabel('Importance')
	xlabel('CF (kHz)')
	set(gca, 'fontsize', fontsize)
	grid on
	%set(gca, 'yscale', 'log')
	if iinstru == 1
	legend('BE', 'BS', 'H', 'F', 'fontsize', legsize)
	end
	clear isMTF
	set(gca, 'FontSize', fontsize)

	%% D. Plot subset
	mean_acc = mean(accuracies,2);
	std_acc = std(accuracies, [], 2);
	nmodels = length(num_neurons);

	nexttile
	errorbar(num_neurons(1:nmodels/2), mean_acc(1:nmodels/2), std_acc(1:nmodels/2)/sqrt(10));
	hold on
	errorbar(num_neurons(nmodels/2+1:end), mean_acc(nmodels/2+1:end), std_acc(nmodels/2+1:end)/sqrt(10))
	xlabel('# Neurons')
	ylabel('Accuracy')
	if iinstru == 1
	hleg = legend('Best', 'Worst', 'fontsize', legsize);
	hleg.ItemTokenSize = [8, 8];
	end
	set(gca, 'FontSize', fontsize)
	box off
	grid on

	%% E. Plot linear

	nexttile
	scatter(10.^(T.response), 10.^(pred_F0), scattersize, 'filled', 'MarkerFaceAlpha',0.3)
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
	plot(x, y, 'r')

	title(['R^2 = ' num2str(r2)])
	xlabel('Actual F0s (Hz)')
	ylabel('Predicted F0s (Hz)')
	set(gca, 'FontSize', fontsize)
	hleg = legend('Fits', 'Unity', 'Fit Line');
	hleg.ItemTokenSize = [8, 8];
	xlim([10^F0s(1) 10^F0s(end)])

	closest_cat = [];
	for ii = 1:length(pred_F0)
		differences = abs(F0s - pred_F0(ii));
		[~, closest_column_index] = min(differences);
		closest_cat(ii) = F0s(closest_column_index);
	end

	% Plot confusion matrix
	nexttile
	imagesc(C)
	title(sprintf('Accuracy = %0.2f%%', accuracy*100))
	set(gca, 'FontSize', fontsize)

%% Save data 

base = getPathsNT;

end

%% Arrange plots

%% Annotate

annotation("textbox", [0.03517 0.602 0.1386 0.1088], "String", "Bassoon",...
	"FontSize", 14, "FontWeight", "bold", "EdgeColor", "none", "Rotation",90)
annotation("textbox", [0.0148 0.1858 0.09778 0.05666], "String", "Oboe",...
	"FontSize", 14, "FontWeight", "bold", "EdgeColor", "none", "Rotation",90)


%% Save figure
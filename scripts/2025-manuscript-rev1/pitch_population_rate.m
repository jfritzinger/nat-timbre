%% plot_pitch_population_rate
clear
save_fig = 1;

%% Set up figure
[base, ~, ~, ppi] = getPathsNT();

figure('Position',[50,50,5.6*ppi,5*ppi])
%tiledlayout(2, 6, 'TileSpacing','tight', 'Padding','compact')
linewidth = 1;
fontsize = 8;
legsize = 7;
titlesize = 10;
labelsize = 12;
scattersize = 5;
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
	load(fullfile(base, 'model_comparisons', ['Pop_Rate_F0_Linear_' target '2.mat']),...
		"pred_F0_test", "T", "r2", "Mdl")
	load(fullfile(base, 'model_comparisons', ['Pop_Rate_F0_' target '_CF2.mat']), ...
		"accuracy", "num_subset")
	accuracies_CF = accuracy;
	load(fullfile(base, 'model_comparisons', ['Pop_Rate_F0_' target '_MTF2.mat']), ...
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

	h(ind(iinstru)+1) = subplot(4, 5, ind(iinstru)+1);
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
	if iinstru == 1 || iinstru == 2
		fprintf('%s SVM accuracy = %0.02f\n', target, pop_rate_F0.validationAccuracy)
	else
		fprintf('%s SVM accuracy = %0.02f\n', target, pop_rate_F0.accuracy)
	end

	%% B. Plot linear

	h(ind(iinstru)+2) = subplot(4, 5, ind(iinstru)+2);
	[T_train, T_test] = splitData_Reps(T);
	scatter(10.^(T_test.Response), 10.^(pred_F0_test), scattersize, 'filled', ...
		'MarkerFaceAlpha',0.5, 'MarkerFaceColor','k')
	set(gca, 'xscale', 'log', 'yscale', 'log')
	hold on
	plot(10.^T_test.Response, 10.^T_test.Response, 'k')

	mdl = fitlm(T_test.Response, pred_F0_test);


	x = linspace(58, 1661, 20);
	p(1) = mdl.Coefficients.Estimate(2,1);
	p(2) = mdl.Coefficients.Estimate(1,1);
	p(3) = mdl.Coefficients.pValue(2);
	y = 10^p(2) .* x.^p(1);
	%y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
	plot(x, y, 'r')
	%plot([10 600], [10 600], 'g', 'LineWidth',1)

	% One semitone above 
	F0s1 = F0s(5:end);
	plot(F0s(1:end-4), F0s1, 'g')
	plot(F0s1, F0s(1:end-4), 'g')


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
	fprintf('%s LR R^2 = %0.02f\n', target, r2)


	%% C. Plot subset CF
	mean_acc = accuracies_CF;
	nmodels = size(accuracies_CF, 2);

	h(ind(iinstru)+3) = subplot(4, 5, ind(iinstru)+3);
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

	% %% D. Plot subset MTF
	% mean_acc = accuracies;
	% nmodels = size(accuracies, 2);
	% 
	% h(ind(iinstru)+4) = subplot(4, 5, ind(iinstru)+4);
	% plot(num_subset(1:nmodels/2), mean_acc(:, 1:nmodels/2));
	% hold on
	% %errorbar(num_neurons(nmodels/2+1:end), mean_acc(nmodels/2+1:end), std_acc(nmodels/2+1:end)/sqrt(10))
	% if iinstru == 1
	% 	hleg = legend(['Best', MTF_names(2:end)], 'fontsize', ...
	% 		legsize, 'box', 'off');
	% 	hleg.ItemTokenSize = [8, 8];
	% elseif iinstru == 3
	% 	xlabel('# Neurons')
	% end
	% set(gca, 'FontSize', fontsize)
	% box off
	% grid on
	% 
	% 	F0s = getF0s(target);
	% F0s = log10(F0s);
	% [sesh, num_data] = getF0Sessions(nat_data, target);
	% T = getF0PopTable(nat_data, target, sesh, F0s, num_data, 'linear', 'Rate');
	% CFs = [nat_data(sesh).CF];
	% MTFs = {nat_data(sesh).MTF};
	% putative = {nat_data(sesh).putative};
	% 
	% avgBeta = Mdl.Beta;  % Average across folds
	% h(ind(iinstru)+3) = subplot(3, 5, ind(iinstru)+3);
	% mdl = fitlm(CFs, avgBeta);
	% x = linspace(300, 14000, 50);
	% y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
	% hleg = legend('', ...
	% 	sprintf('p=%0.04f',mdl.Coefficients{2,4}), 'fontsize', legsize, ...
	% 	'location', 'northwest', 'box', 'off');

end

%% 

load(fullfile(base, 'model_comparisons', 'Pop_Rate_F0_Linear_F0s.mat'),...
	"results")
r2_bassoon = [results(1,:).r2];
r2_oboe = [results(2,:).r2];
r2_invariant = [results(3,:).r2];
mse_bassoon = [results(1,:).test_mse];
mse_oboe = [results(2,:).test_mse];
mse_invariant = [results(3,:).test_mse];
conditions = {'Bassoon', 'Oboe', 'Both'};

h(16) = subplot(4, 5, 16);
hold on
swarmchart(ones(500,1), r2_bassoon, scattersize,'filled', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha',0.5)
swarmchart(ones(500,1)*2, r2_oboe, scattersize,'filled', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha',0.5)
swarmchart(ones(500,1)*3, r2_invariant, scattersize,'filled', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha',0.5)
set(gca, 'FontSize', fontsize)
box off
grid on
ylim([0 1])
xticks(1:3)
ylabel('R^2')
xticklabels(conditions)

h(17) = subplot(4, 5, 17);
hold on
swarmchart(ones(500,1), mse_bassoon, scattersize,'filled', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha',0.5)
swarmchart(ones(500,1)*2, mse_oboe, scattersize,'filled', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha',0.5)
swarmchart(ones(500,1)*3, mse_invariant, scattersize,'filled', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha',0.5)
set(gca, 'FontSize', fontsize)
box off
grid on
xticks(1:3)
ylabel('MSE')
xticklabels(conditions)


%% 

load(fullfile(base, 'model_comparisons', 'Pop_Rate_F0_Linear_Instruments.mat'),...
	"results")
r2_bassoon_train = [results(1,:).r2];
r2_oboe_train = [results(2,:).r2];
mse_bassoon_train = [results(1,:).test_mse];
mse_oboe_train = [results(2,:).test_mse];
conditions = {'Oboe', 'Bassoon'};

h(18) = subplot(4, 5, 18);
hold on
swarmchart(ones(500,1), r2_bassoon_train, scattersize, 'filled', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha',0.5)
swarmchart(ones(500,1)*2, r2_oboe_train, scattersize,'filled', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha',0.5)
set(gca, 'FontSize', fontsize)
ylim([0 1])
box off
grid on
xticks(1:2)
ylabel('R^2')
xticklabels(conditions)
xlabel('Test Instrument')

h(19) = subplot(4, 5, 19);
hold on
swarmchart(ones(500,1), mse_bassoon_train, scattersize,'filled', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha',0.5)
swarmchart(ones(500,1)*2, mse_oboe_train, scattersize,'filled', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha',0.5)
set(gca, 'FontSize', fontsize)
box off
grid on
xticks(1:2)
ylabel('MSE')
xticklabels(conditions)
xlabel('Test Instrument')

%% Arrange plots

left = [0.100    0.42    0.64    0.83];
bottom = [fliplr(linspace(0.38, 0.8, 3)) 0.09];
width = 0.15;
height = 0.15;

for ii = 1:3
	set(h(ind(ii)+1), 'position', [left(1) bottom(ii) width height]);
	set(h(ind(ii)+2), 'position', [left(2) bottom(ii) width height]);
	set(h(ind(ii)+3), 'position', [left(3) bottom(ii) width height]);
	set(h(ind(ii)+4), 'position', [left(4) bottom(ii) width height]);
	%set(h(ind(ii)+5), 'position', [left(5) bottom(ii) width height]);
end

left = linspace(0.08, 0.8, 4);
width = 0.18;
height = 0.15;

set(h(16), 'position', [left(1) bottom(4) width height]);
set(h(17), 'position', [left(2) bottom(4) width height]);
set(h(18), 'position', [left(3) bottom(4) width height]);
set(h(19), 'position', [left(4) bottom(4) width height]);

%% Annotate

annotation("textbox", [0.09 0.78 0.1386 0.1088], "String", "Bassoon",...
	"FontSize", titlesize, "FontWeight", "bold", "EdgeColor", "none", "Rotation",90)
annotation("textbox", [0.04 0.6 0.09778 0.05666], "String", "Oboe",...
	"FontSize", titlesize, "FontWeight", "bold", "EdgeColor", "none", "Rotation",90)
annotation("textbox", [0.09 0.4 0.1386 0.1088], "String", "Both",...
	"FontSize", titlesize, "FontWeight", "bold", "EdgeColor", "none", "Rotation",90)

annotation("textbox", [0.17 0.2 0.50 0.1088], "String", "Test on held-out F0s",...
	"FontSize", titlesize, "FontWeight", "bold", "EdgeColor", "none")
annotation("textbox", [0.61 0.2 0.50 0.1088], "String", "Test on opposite instrument",...
	"FontSize", titlesize, "FontWeight", "bold", "EdgeColor", "none")


left = [0.0900    0.42    0.64    0.83];
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
annotation('textbox',[0.04 0.25 0.071 0.058],...
	'String','E','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[0.53 0.25 0.071 0.058],...
	'String','F','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
%% Save figure

if save_fig == 1
	filename = 'fig8_pitch_pop_rate';
	save_figure(filename)
end

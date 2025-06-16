%% timbre_pop_betas 
clear 
save_fig = 1;

%% Load in data 

[base, datapath, savepath, ppi] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')
load(fullfile(base, 'model_comparisons', 'Pop_Rate_Timbre_All.mat'), ...
	"pop_rate_timbre")

sesh = [];
for ii = 1:length(nat_data)
	rate = nat_data(ii).bass_rate;
	rate2 = nat_data(ii).oboe_rate;
	if ~isempty(rate) && ~isempty(rate2)
		sesh = [sesh ii];
	end
end
num_data = length(sesh);

beta_weights = pop_rate_timbre.trainedClassifier.ClassificationSVM.Beta;

%% Create figure 

figure('Position',[50, 50, 6.8*ppi, 4*ppi])
tiledlayout(1, 2, 'Padding','compact')
linewidth = 2;
fontsize = 18;
titlesize = 20;
legsize = 20;
scattersize = 40;
colorsData = {"#0072BD", "#D95319"};
colorsMTF = {'#648FFF', '#DC267F', '#785EF0', '#FFB000'};
colorsTimbre = '#1b9e77';

%% E. Beta weights, CF group

nexttile

CFs = pop_rate_timbre.CFs;
scatter(CFs, beta_weights, scattersize, 'filled', 'MarkerEdgeColor','k', ...
	'MarkerFaceColor',colorsTimbre, 'MarkerFaceAlpha',0.5);
hold on
set(gca, 'xscale', 'log')

mdl = fitlm(CFs, beta_weights);
x = linspace(300, 14000, 50);
y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
plot(x, y, '--k', 'linewidth', 2)
hleg = legend('', ...
	sprintf('p=%0.04f',mdl.Coefficients{2,4}), 'fontsize', legsize, ...
	'location', 'northwest', 'box', 'off');
hleg.ItemTokenSize = [8, 8];
ylabel('Beta Weights')
xlabel('CFs')
ylim([-1.8 3])
xticks([200 500 1000 2000 5000 10000])
xticklabels([0.2 0.5 1 2 5 10])
xtickangle(0)
set(gca, 'fontsize', fontsize)
grid on

%% D. Beta weights, MTF types 

nexttile
CFs = pop_rate_timbre.CFs;
MTFs = pop_rate_timbre.MTF;
MTF_types = unique(MTFs);

[weights_ordered, order_ind] = sort(beta_weights);
hold on
for iMTF = 1:4
	if iMTF == 4
		ind = strcmp(MTFs, MTF_types{iMTF}) | strcmp(MTFs, MTF_types{iMTF+1});
	else
		ind = strcmp(MTFs, MTF_types{iMTF});
	end
	[weights_ordered, order_ind] = sort(beta_weights(ind));
	num_units = length(weights_ordered);

	RGB = hex2rgb(colorsMTF{iMTF});
	swarmchart(ones(num_units, 1)*iMTF, weights_ordered, scattersize, ...
		RGB, "filled", 'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.5)
	ylim([-2 2])

	mean_vals(iMTF) = mean(weights_ordered);
	std_vals(iMTF) = std(weights_ordered)/sqrt(length(weights_ordered));
end
errorbar(1:4, mean_vals, std_vals, 'k')
xticks(1:4)
xticklabels({'BE', 'BS', 'F', 'H'})
xlabel('MTF Groups')
ylabel('Beta Weights')

tableMTF = table(MTFs', beta_weights);
% anova(tableMTF, 'beta_weights')
% [~,~,stats] = anova1(beta_weights, MTFs);
% [c,~,~,gnames] = multcompare(stats);
set(gca, 'fontsize', fontsize)
grid on

%% Save figure 

if save_fig == 1
	filename = 'timbre_pop_betas';
	save_figure_MARC(filename)
end

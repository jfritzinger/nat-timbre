%% plot_pca_vs_CF
clear

%% Load in data_NT_matrix

% Load in spreadsheet
base = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre';
filepath = 'data/model_comparisons';
filename = 'Data_NT_Matrix.mat';
load(fullfile(base,filepath,filename), "data_bassoon", "data_oboe");

%% Plot

figure('Position',[1771,586,509,712])
tiledlayout(2, 1)
textsize = 14;

% Bassoon
rate_matrix = data_bassoon.rates_z;
CFs = data_bassoon.CFs;
nexttile
[~, score] = pca(rate_matrix);
scatter(CFs, score(:,1), 'filled', 'MarkerEdgeColor','k')
title('Bassoon, PCA 1 vs CF')
xlabel('CF')
set(gca, 'xscale', 'log')
r = corrcoef(log(CFs), score(:,1));
r2 = r(1, 2)^2;
text(0.05, 0.95, ['R^2=' num2str(round(r2, 2))], 'Units', 'normalized', ...
	'VerticalAlignment', 'top', 'FontSize',16)
set(gca, 'fontsize', textsize)
ylabel('Component 1')

% Oboe
rate_matrix = data_oboe.rates_z;
CFs = data_oboe.CFs;
nexttile
[~, score] = pca(rate_matrix);
scatter(CFs, score(:,1), 'filled', 'MarkerEdgeColor','k')
title('Oboe, PCA 1 vs CF')
xlabel('CF')
set(gca, 'xscale', 'log')
r = corrcoef(log(CFs), score(:,1));
r2 = r(1, 2)^2;
text(0.05, 0.95, ['R^2=' num2str(round(r2, 2))], 'Units', 'normalized', ...
	'VerticalAlignment', 'top', 'FontSize',16)
set(gca, 'fontsize', textsize)
ylabel('Component 1')

%% Split bassoon into voice pitch and energy regions

figure('Position',[1771,586,509,712])
tiledlayout(2, 1)
textsize = 14;

% Bassoon
% 27:end is 277 Hz to 587 Hz
% 3:26 is 64 Hz to 262 Hz
for ii = 1:2
	if ii == 1
		rate_matrix = data_bassoon.rates_z(:,3:26);
	else
		rate_matrix = data_bassoon.rates_z(:,27:end);
	end
	CFs = data_bassoon.CFs;
	nexttile
	[~, score] = pca(rate_matrix);
	scatter(CFs, score(:,1), 'filled', 'MarkerEdgeColor','k')
	xlabel('CF')
	set(gca, 'xscale', 'log')
	r = corrcoef(log(CFs), score(:,1));
	r2 = r(1, 2)^2;
	text(0.05, 0.95, ['R^2=' num2str(round(r2, 2))], 'Units', 'normalized', ...
		'VerticalAlignment', 'top', 'FontSize',16)
	set(gca, 'fontsize', textsize)
	ylabel('Component 1')
	if ii == 1
		title('Bassoon, Voice pitch (64-262Hz)')
	else
		title('Bassoon, High pitch (277-587Hz)')
	end
end

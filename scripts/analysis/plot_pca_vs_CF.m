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

%% PCA VARIANCE OF TUNING ANALYSIS

rate_matrix = data_bassoon.rates_z;
allTuningMats = rate_matrix';
[height, num_mats] = size(allTuningMats); %measure sizes of image stack
vectMats = reshape(allTuningMats, height, num_mats)';  % vectorize each tuning matrix
X_centered = vectMats - mean(vectMats, 1); %subtract off means
[coeff, score, latent] = pca(X_centered); %run PCA

% Plot variance explained
figure;
explained = 100*latent/sum(latent);
plot(cumsum(explained),'k-'); hold on;
plot(cumsum(explained),'k.'); hold off;
xlabel('Number of Components');
ylabel('Variance Explained (%)');
title('Cumulative Variance Explained');
grid on

% Visualize top 8 components
num_components_to_show = 16; 
figure; tiledlayout(4,4, 'TileSpacing','compact', 'Padding','compact');
for n = 1:num_components_to_show
    nexttile(n);
    pcMat = coeff(:,n);
    plot(pcMat); %clim([-0.5 0.8]);
    title(['PC ' num2str(n)]);
end

% Find images that are most aligned with each principal component
figure; tiledlayout(4,4, 'TileSpacing','compact', 'Padding','compact');
for n = 1:num_components_to_show
    [~, max_idx] = max(abs(score(:,n)));
    nexttile(n)
    orig_image = allTuningMats(:,max_idx);
    plot(orig_image); %clim([-4 8]);
    title(['Neuron most aligned with PC ' num2str(n)],'FontSize',8);
end
%% timbre_pop_timing
clear 
save_fig = 1;

%% Load in data 

[base, ~, ~, ppi] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'pop_timing_timbre.mat'), ...
	"accur_all","C_all", "num_neurons", "mean_acc", "std_acc")


%% Plot 

figure('Position',[50 50, 5*ppi, 3*ppi])
linewidth = 2;
fontsize = 18;
titlesize = 20;
legsize = 20;
scattersize = 40;

nmodels = length(num_neurons);
nexttile
errorbar(num_neurons(1:nmodels/2), mean_acc(1:nmodels/2),...
	std_acc(1:nmodels/2)/sqrt(1), 'Color','#1b9e77', 'LineWidth',2);
hold on
% errorbar(num_neurons(nmodels/2+1:end), mean_acc(nmodels/2+1:end), ...
% 	std_acc(nmodels/2+1:end)/sqrt(10), 'Color','k', 'LineWidth',2)
xlabel('# Neurons in Model')
ylabel('Accuracy')
box off
grid on
set(gca, 'fontsize', fontsize)
yticks(0.8:0.05:1)
xticks(0:10:40)

% Save figure 
if save_fig == 1
	filename = 'timbre_pop_timing';
	save_figure_MARC(filename)
end

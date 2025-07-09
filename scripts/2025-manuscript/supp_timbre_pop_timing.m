%% supp_timbre_pop_timing

%% Load in data 

[base, ~, ~, ppi] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'pop_timing_timbre.mat'), ...
	"accur_all","C_all", "num_neurons", "mean_acc", "std_acc")
load(fullfile(base,'model_comparisons', 'Data_NT_3.mat'), 'nat_data')


%% Plot 

figure('Position',[50,50,2.5*ppi,2*ppi])
linewidth = 1;
fontsize = 8;
legsize = 7;
titlesize = 10;
labelsize = 12;
scattersize = 10;

nmodels = length(num_neurons);
nexttile
errorbar(num_neurons(1:nmodels/2), mean_acc(1:nmodels/2), ...
	std_acc(1:nmodels/2)/sqrt(5), 'Color','#1b9e77', ...
	'LineWidth',linewidth);
hold on
errorbar(num_neurons(nmodels/2+1:end), mean_acc(nmodels/2+1:end),...
	std_acc(nmodels/2+1:end)/sqrt(10), 'Color','k', 'LineWidth',linewidth)
xlabel('Number of Neurons in Model')
ylabel('Accuracy')
legend('Best', 'Worst')
box off
grid on
set(gca, 'fontsize', fontsize)

%% Save figure

save_fig = 1;
if save_fig == 1
	filename = 'supp1_timbrepoptiming';
	save_figure(filename)
end
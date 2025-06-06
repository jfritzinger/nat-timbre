%% plot_neuron_all
clear 

%% Load in rate 


[base, datapath, savepath, ppi] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Neuron_Rate_All.mat'), ...
	"neuron_rate_all")
load(fullfile(base, 'model_comparisons', 'Neuron_Time_All.mat'), ...
	"neuron_time_all")

%% 

accuracy_rate = [neuron_rate_all.accuracy];
accuracy_time = [neuron_time_all.accuracy];

% Plot rate 
figure
tiledlayout(1, 3)
nexttile
edges = linspace(0, 1, 31);
histogram(accuracy_rate, edges)
hold on
xline(1/75)
xlabel('Accuracy')
title('Rate Model')
ylabel('# Neurons')

% Plot timing 
nexttile
edges = linspace(0, 1, 31);
histogram(accuracy_time, edges);
hold on
xline(1/75)
xlabel('Accuracy')
title('PSTH Model')
ylabel('# Neurons')

% Plot rate against timing 
nexttile
scatter(accuracy_rate, accuracy_time, 10, 'filled', 'MarkerEdgeColor','k')
xlabel('Rate Accuracy')
ylabel('Time Accuracy')
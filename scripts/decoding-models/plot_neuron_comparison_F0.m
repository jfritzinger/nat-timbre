%% plot_neuron_comparison_F0.m
clear

%% Load in data 

target = 'Oboe';

[base, datapath, savepath, ppi] = getPathsNT();
load(fullfile(base, 'model_comparisons', ['Neuron_Time_F0_' target '.mat']), ...
	"neuron_time_F0")
load(fullfile(base, 'model_comparisons', ['Neuron_Rate_F0_' target '.mat']), ...
	"neuron_rate_F0")

%% Plot accuracy of each neuron
figure('Position',[560,618,798,230])
tiledlayout(1, 3)
linewidth = 1;

accuracy_time = [neuron_time_F0.accuracy]*100;
accuracy_rate = [neuron_rate_F0.accuracy]*100;

% Plot overall histogram
nexttile
hold on
edges = linspace(0, 100, 30);
histogram(accuracy_rate,edges, 'FaceColor','r')
histogram(accuracy_time,edges, 'FaceColor','b')
chance = 1/size(neuron_time_F0(1).C, 1)*100;
xline(chance, 'k', 'LineWidth',linewidth)
xline(mean(accuracy_time), 'b', 'LineWidth',2)
xline(mean(accuracy_rate), 'r', 'LineWidth',2)
ylabel('# Neurons')
xlabel('Prediction Accuracy (%)')
title([target ' Prediction of F0'])
xlim([0 65])
legend('', '', sprintf('Chance = %.2f%%', chance), ...
	sprintf('Mean = %.2f%%', mean(accuracy_time)), ...
	sprintf('Mean = %.2f%%', mean(accuracy_rate)))

% Best neuron confusion matrices 

% Find indices of the best neurons 
[temp,originalpos] = sort(accuracy_time, 'descend' );
n = temp(1:3);
best_ind=originalpos(1:3);

% Plot confusion matrix of best neuron
putatives = neuron_time_F0(best_ind(1)).putative;
nexttile
confusionchart(neuron_time_F0(best_ind(1)).C)
title(sprintf('Timing, %s, %0.02f%%', putatives, temp(1)))

nexttile
confusionchart(neuron_rate_F0(best_ind(1)).C)
title(sprintf('Rate, %s, %0.02f%%', putatives, accuracy_rate(best_ind(1))))

%% plot_timbre_comparison
clear

%% Get best single-unit timbre neurons 

[base, ~, ~, ppi] = getPathsNT();
filepath = fullfile(base, 'model_comparisons','Neuron_Rate_Timbre_All.mat');
load(filepath, "neuron_rate_timbre")
accuracy_rate = [neuron_rate_timbre.accuracy]*100;


%% Get best population timbre neurons 

load(fullfile(base, 'model_comparisons', 'Pop_Rate_Timbre_All.mat'), ...
	"pop_rate_timbre")
beta_weights = pop_rate_timbre.trainedClassifier.ClassificationSVM.Beta;
CFs = pop_rate_timbre.CFs;

%% Are they the same neurons? 

figure
hold on
scatter(accuracy_rate, abs(beta_weights), 'filled', 'MarkerEdgeColor','k')
xlabel('Accuracy for Single-Unit Rate')
xlim([40 100])
xticks(0:5:100)

ylabel('|Beta Weights| for Population Rate')
ylim([0 2])
yticks(0:0.2:5)

mdl = fitlm(accuracy_rate, abs(beta_weights));
x = linspace(0, 100, 20);
y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
plot(x, y, 'r')
hleg = legend('Neuron', 'p=7.1671e-25');
hleg.ItemTokenSize = [8, 8];
title('Single Neuron vs Population Betas')
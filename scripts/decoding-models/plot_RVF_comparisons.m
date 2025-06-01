
clear 

%% Load in timbre single-neuron results 

[base, datapath, savepath, ppi] = getPathsNT();
filepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/model_comparisons';
load(fullfile(filepath, 'Data_NT_3.mat'), 'nat_data')

% Load in data 
filename = 'Neuron_Rate_Timbre_All.mat';
filepath = fullfile(base, 'model_comparisons',filename);
load(filepath, "neuron_rate_timbre")


%% Get subset of units 

sesh = [];
for ii = 1:length(nat_data)
	rate = nat_data(ii).bass_rate;
	rate2 = nat_data(ii).oboe_rate;
	if ~isempty(rate) && ~isempty(rate2)
		sesh = [sesh ii];
	end
end
num_data = length(sesh);

PCA_scores = [nat_data.RVF_PC2];
PC2_score = PCA_scores(sesh)';
accuracy = [neuron_rate_timbre.accuracy];

%%

figure
tiledlayout(1, 3)

nexttile
scatter(PC2_score, accuracy, 'filled')
xlabel('PC2 RVF score')
ylabel('Accuracy')
hold on

mdl = fitlm(PC2_score, accuracy);
x = linspace(-2, 1, 20);
y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
plot(x, y, 'r')
hleg = legend('Neuron', ...
	sprintf('y = %0.2f*x+%0.2f, p=%0.04f', mdl.Coefficients{2, 1}, ...
	mdl.Coefficients{1, 1},mdl.Coefficients{2,4}));
title('Neuron Rate')

%% Timing 

% Load in data 
filename = 'Neuron_Time_Timbre_All.mat';
filepath = fullfile(base, 'model_comparisons',filename);
load(filepath, "neuron_time_timbre")
accuracy = [neuron_time_timbre.accuracy];

nexttile
scatter(PC2_score, accuracy, 'filled')
xlabel('PC2 RVF score')
ylabel('Accuracy')
hold on

mdl = fitlm(PC2_score, accuracy);
x = linspace(-2, 1, 20);
y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
plot(x, y, 'r')
hleg = legend('Neuron', ...
	sprintf('y = %0.2f*x+%0.2f, p=%0.04f', mdl.Coefficients{2, 1}, ...
	mdl.Coefficients{1, 1},mdl.Coefficients{2,4}));
title('Neuron Timing')

%% Load in timbre population results 

load(fullfile(base, 'model_comparisons', 'Pop_Rate_Timbre_All.mat'), ...
	"pop_rate_timbre")

% Plot PCA2 score vs beta weights 
beta = pop_rate_timbre.trainedClassifier.ClassificationSVM.Beta;

nexttile
scatter(PC2_score, beta, 'filled')
xlabel('PC2 RVF score')
ylabel('Beta')
hold on

mdl = fitlm(PC2_score, beta);
x = linspace(-2, 1, 20);
y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
plot(x, y, 'r')
hleg = legend('Neuron', ...
	sprintf('y = %0.2f*x+%0.2f, p=%0.04f', mdl.Coefficients{2, 1}, ...
	mdl.Coefficients{1, 1},mdl.Coefficients{2,4}));
title('Population')

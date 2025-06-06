%% timbre_population
clear 
save_fig = 0;

%% Load in data 

[base, datapath, savepath, ppi] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')
load(fullfile(base, 'model_comparisons', 'Pop_Rate_Timbre_All.mat'), ...
	"pop_rate_timbre")
load(fullfile(base, 'model_comparisons','Neuron_Rate_Timbre_All.mat'),...
	"neuron_rate_timbre")
accuracy_rate = [neuron_rate_timbre.accuracy]*100;

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
%beta_weights = pop_rate_timbre.imp.ImportanceMean;

%% Create figure 

figure('Position',[50, 50, 6*ppi, 5*ppi])
scattersize = 5;
titlesize = 9;
fontsize = 8;
labelsize = 12;
legsize = 7;


%% A. Example 

h = plot_example_populations_timbre();

%% B. Accuracies for shuffled data vs regular

h(4) = subplot(4, 3, 5);
edges = linspace(0, 1, 41);
hold on
histogram(pop_rate_timbre.shuffled_accuracy, edges)
xline(pop_rate_timbre.accuracy, '--r')
xline(0.5, 'k')
xlim([0 1])
xlabel('Accuracy (%)')
ylabel('# Trials')
hleg = legend(['Shuffled' newline 'Data'], ['Model' newline 'Accuracy']...
	, 'Chance', 'fontsize', legsize, 'location', 'best', 'box', 'off');
hleg.ItemTokenSize = [8, 8];
set(gca, 'fontsize', fontsize)
grid on

%% C. Beta weights vs single neuron accuracy 
h(5) = subplot(4, 3, 6);

hold on
scatter(accuracy_rate, abs(beta_weights),scattersize, 'filled', 'MarkerEdgeColor','k')
xlabel('Single-Unit Rate Accuracy')
xlim([40 100])
xticks(0:20:100)

ylabel('|Beta Weights|')
ylim([0 2.5])
yticks(0:0.5:5)

mdl = fitlm(accuracy_rate, abs(beta_weights));
x = linspace(0, 100, 20);
y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
plot(x, y, 'r')
hleg = legend('Neuron',sprintf('p=%0.04f',mdl.Coefficients{2,4}), ...
	'fontsize', legsize, 'location', 'northwest', 'box', 'off');
hleg.ItemTokenSize = [8, 8];
set(gca, 'fontsize', fontsize)
grid on


%% E. Beta weights, CF group

h(6) = subplot(4, 3, 7);
% CF_groups = [0, 2000; 2000, 4000; 4000, 14000];
% CF_names = {'Low', 'Medium', 'High'};
% hold on
% mean_vals = zeros(1, 3);
% std_vals = zeros(1,3);
% for iCF = 1:3
% 	ind = CFs > CF_groups(iCF, 1) & CFs < CF_groups(iCF, 2);
% 	[weights_ordered, order_ind] = sort(beta_weights(ind));
% 	num_units = length(weights_ordered);
% 
% 	swarmchart(ones(num_units, 1)*iCF, weights_ordered, scattersize)
% 	ylim([-2 2])
% 
% 	mean_vals(iCF) = mean(weights_ordered);
% 	std_vals(iCF) = std(weights_ordered)/sqrt(length(weights_ordered));
% end
% errorbar(1:3, mean_vals, std_vals, 'k')
% xticks(1:3)
% xticklabels(CF_names)
% ylabel('Beta Weights')
% 
% tableMTF = table(CFs', beta_weights);
% anova(tableMTF, 'beta_weights')
% % [~,~,stats] = anova1(beta_weights, CFs);
% % [c,~,~,gnames] = multcompare(stats);
% set(gca, 'fontsize', fontsize)

CFs = pop_rate_timbre.CFs;
scatter(CFs, beta_weights, scattersize, 'filled', 'MarkerEdgeColor','k');
hold on
set(gca, 'xscale', 'log')

mdl = fitlm(CFs, beta_weights);
x = linspace(300, 14000, 50);
y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
plot(x, y, 'r')
hleg = legend('Neuron', ...
	sprintf('p=%0.04f',mdl.Coefficients{2,4}), 'fontsize', legsize, ...
	'location', 'northwest', 'box', 'off');
hleg.ItemTokenSize = [8, 8];
ylabel('Beta Weights')
xlabel('CFs')
ylim([-1.8 3])
xticks([200 500 1000 2000 5000 10000])
xticklabels([0.2 0.5 1 2 5 10])
set(gca, 'fontsize', fontsize)
grid on

%% D. Beta weights, MTF types 

h(7) = subplot(4, 3, 8);

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

	swarmchart(ones(num_units, 1)*iMTF, weights_ordered, scattersize)
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
anova(tableMTF, 'beta_weights')
% [~,~,stats] = anova1(beta_weights, MTFs);
% [c,~,~,gnames] = multcompare(stats);
set(gca, 'fontsize', fontsize)
grid on

%% F. Beta weights vs PCA for RVF 

h(8) = subplot(4, 3, 9);

PCA_scores = [nat_data.RVF_PC2];
PC2_score = PCA_scores(sesh)';

% Plot PCA2 score vs beta weights 
scatter(PC2_score, beta_weights, scattersize, 'filled', 'MarkerEdgeColor','k')
xlabel('PC2 RVF score')
ylabel('Beta Weights')
hold on

mdl = fitlm(PC2_score, beta_weights);
x = linspace(-2, 1, 20);
y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
plot(x, y, 'r')
hleg = legend('Neuron', ...
	sprintf('p=%0.04f',mdl.Coefficients{2,4}), 'fontsize', legsize, ...
	'location', 'northwest', 'box', 'off');
hleg.ItemTokenSize = [8, 8];
set(gca, 'fontsize', fontsize)
ylim([-1.8 3])
grid on

%% H. Subset CF groups 
h(9) = subplot(4, 3, 10);

load(fullfile(base, 'model_comparisons', 'Pop_Rate_Timbre_Subset_CF.mat'), ...
	"accuracies", "CF_names")

hold on
for iCF = 1:7
	swarmchart(ones(size(accuracies, 2), 1)*iCF, accuracies(iCF,:), scattersize-4)
end
xlabel('CF Groups')
xticks(1:7)
xticklabels(CF_names)
max_acc = max(accuracies, [], 2);
plot(1:7, max_acc, 'k')
mean_acc = median(accuracies, 2);
plot(1:7, mean_acc, 'k')
set(gca, 'fontsize', fontsize)
ylabel('Model Accuracy')
grid on

%% I. Subset MTF groups 
h(10) = subplot(4, 3, 11);

load(fullfile(base, 'model_comparisons', 'Pop_Rate_Timbre_Subset_MTF.mat'), ...
	"accuracy_all", "MTF_names")

hold on
for iCF = 1:5
	swarmchart(ones(size(accuracy_all, 2), 1)*iCF, accuracy_all(iCF,:), scattersize-4)
end
xlabel('MTF Groups')
xticks(1:5)
xticklabels(MTF_names)
ylim([0.78 1])
max_acc = max(accuracy_all, [], 2);
plot(1:5, max_acc, 'k')
mean_acc = median(accuracy_all, 2);
plot(1:5, mean_acc, 'k')
set(gca, 'fontsize', fontsize)
ylabel('Model Accuracy')
grid on

%% G. Subset accuracies 
h(11) = subplot(4, 3, 12);

load(fullfile(base, 'model_comparisons', 'Pop_Rate_Timbre_Best.mat'), ...
	"pop_rate_timbre_best")

nneurons = [1:4 5:5:50];
accuracy_bad = [pop_rate_timbre_best(15:28).accuracy];
plot(nneurons, accuracy_bad);

hold on 
accuracy_good = [pop_rate_timbre_best(1:14).accuracy];
plot(nneurons, accuracy_good);
xlabel('# Neurons in Model')
ylabel('Model Accuracy')
grid on
box off
hleg = legend('Worst', 'Best', 'fontsize', legsize, 'location', 'best', 'box', 'off');
hleg.ItemTokenSize = [8, 8];
set(gca, 'fontsize', fontsize)

%% Set position

left = linspace(0.1, 0.74, 3);
bottom = linspace(0.13, 0.75, 3);
height = 0.21;
width = 0.23;

set(h(3), 'Position', [left(1) bottom(3) width height/3])
set(h(2), 'Position', [left(1) bottom(3)+height/3 width height/3])
set(h(1), 'Position', [left(1) bottom(3)+height/3*2 width height/3])

set(h(4), 'position', [left(2) bottom(3) width height])
set(h(5), 'position', [left(3) bottom(3) width height])
set(h(6), 'position', [left(1) bottom(2) width height])
set(h(7), 'position', [left(2) bottom(2) width height])
set(h(8), 'position', [left(3) bottom(2) width height])
set(h(9), 'position', [left(1) bottom(1) width height])
set(h(10), 'position', [left(2) bottom(1) width height])
set(h(11), 'position', [left(3) bottom(1) width height])


%% Annotate 

left = linspace(0.01, 0.68, 3);
bottom = linspace(0.36, 0.97, 3);

label = {'A', 'B', 'C'};
for ii = 1:3
	annotation('textbox',[left(ii) bottom(3) 0.0826 0.0385],'String',label{ii},...
		'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
end
label = {'D', 'E', 'F'};
for ii = 1:3
	annotation('textbox',[left(ii) bottom(2) 0.0826 0.0385],'String',label{ii},...
		'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
end
label = {'G', 'H', 'I'};
for ii = 1:3
	annotation('textbox',[left(ii) bottom(1) 0.0826 0.0385],'String',label{ii},...
		'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
end

%% Save figure 

if save_fig == 1
	filename = 'fig3_timbre_population_rate';
	save_figure(filename)
end
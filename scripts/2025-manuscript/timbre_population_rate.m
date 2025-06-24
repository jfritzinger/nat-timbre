%% timbre_population
clear 
save_fig = 1;

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
colorsData = {"#0072BD", "#D95319"};
colorsMTF = {'#648FFF', '#DC267F', '#785EF0', '#FFB000'};
colorsTimbre = '#1b9e77';
linewidth = 1;

%% A. Example 

h = plot_example_populations_timbre();

%% B. Accuracies for shuffled data vs regular

h(4) = subplot(4, 3, 5);
edges = linspace(0, 1, 41);
hold on
histogram(pop_rate_timbre.shuffled_accuracy, edges, ...
	'FaceColor', 'k')
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
scatter(accuracy_rate, abs(beta_weights),scattersize, 'filled', ...
	'MarkerEdgeColor','k', 'MarkerFaceColor',colorsTimbre)
xlabel('Single-Unit Rate Accuracy')
xlim([40 100])
xticks(0:20:100)

ylabel('|Beta Weights|')
ylim([0 2.5])
yticks(0:0.5:5)

mdl = fitlm(accuracy_rate, abs(beta_weights));
x = linspace(0, 100, 20);
y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
plot(x, y, ':k')
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
scatter(CFs, beta_weights, scattersize, 'filled', 'MarkerEdgeColor','k', ...
	'MarkerFaceColor',colorsTimbre);
hold on
set(gca, 'xscale', 'log')

mdl = fitlm(CFs, beta_weights);
x = linspace(300, 14000, 50);
y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
plot(x, y, ':k')
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

	RGB = hex2rgb(colorsMTF{iMTF});
	swarmchart(ones(num_units, 1)*iMTF, weights_ordered, scattersize, RGB)
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
%anova(tableMTF, 'beta_weights')
% [~,~,stats] = anova1(beta_weights, MTFs);
% [c,~,~,gnames] = multcompare(stats);
set(gca, 'fontsize', fontsize)
grid on

%% F. Beta weights vs PCA for RVF 
% 
% h(8) = subplot(4, 3, 9);
% 
% PCA_scores = [nat_data.RVF_PC2];
% PC2_score = PCA_scores(sesh)';
% 
% % Plot PCA2 score vs beta weights 
% scatter(PC2_score, beta_weights, scattersize, 'filled', 'MarkerEdgeColor','k', ...
% 	'MarkerFaceColor',colorsTimbre)
% xlabel('PC2 RVF score')
% ylabel('Beta Weights')
% hold on
% 
% mdl = fitlm(PC2_score, beta_weights);
% x = linspace(-2, 1, 20);
% y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
% plot(x, y, ':k')
% hleg = legend('Neuron', ...
% 	sprintf('p=%0.04f',mdl.Coefficients{2,4}), 'fontsize', legsize, ...
% 	'location', 'northwest', 'box', 'off');
% hleg.ItemTokenSize = [8, 8];
% set(gca, 'fontsize', fontsize)
% ylim([-1.8 3])
% grid on

%% G. Subset CF groups 
h(10) = subplot(4, 3, 11);

colorsCF = [27/256, 158/256, 119/256;0.0660 0.4430 0.7450; 0.8660 0.3290 0.0000;0.9290 0.6940 0.1250];
load(fullfile(base, 'model_comparisons', 'Pop_Rate_Timbre_Subset_MTF2.mat'), ...
	"accuracy", "std_acc", "MTF_names")
acc_MTF = accuracy;
std_MTF = std_acc;
load(fullfile(base, 'model_comparisons', 'Pop_Rate_Timbre_Subset_CF2.mat'), ...
	"accuracy", "std_acc", "CF_groups")
num_subset = [1, 2, 3, 4, 5:5:90];

colorsMTF = {'#1b9e77', '#648FFF', '#DC267F', '#785EF0', '#FFB000'};
hold on
for iMTF = 1:5
	plot(num_subset, acc_MTF(iMTF,:), 'Color', colorsMTF{iMTF}, 'linewidth', linewidth)
end
ylabel('Accuracy')
xlabel('# Neurons in Model')
hleg = legend({'Best Neurons', MTF_names{2:5}}, 'Location','southeast');
hleg.ItemTokenSize = [8, 8];

grid on 
box off
set(gca, 'fontsize', fontsize)
xticks(0:20:100)
xtickangle(0)


% load(fullfile(base, 'model_comparisons', 'Pop_Rate_Timbre_Subset_CF.mat'), ...
% 	"accuracies", "CF_names")
% 
% hold on
% for iCF = 1:7
% 	swarmchart(ones(size(accuracies, 2), 1)*iCF, accuracies(iCF,:), scattersize-4)
% end
% xlabel('CF Groups')
% xticks(1:7)
% xticklabels(CF_names)
% max_acc = max(accuracies, [], 2);
% plot(1:7, max_acc, 'k')
% mean_acc = median(accuracies, 2);
% plot(1:7, mean_acc, 'k')
% set(gca, 'fontsize', fontsize)
% ylabel('Model Accuracy')
% grid on

%% H. Subset accuracies 
h(9) = subplot(4, 3, 10);

hold on
for i = 1:4
	plot(num_subset, accuracy(i,:), 'linewidth', linewidth, 'Color',colorsCF(i,:))
end
ylabel('Accuracy')
xlabel('# Neurons in Model')
hleg = legend('Best Neurons', 'CF = 0-2 kHz', 'CF = 2-4 kHz', ...
	'CF = 4-8 kHz', 'Location','southeast');
hleg.ItemTokenSize = [8, 8];

grid on 
box off
set(gca, 'fontsize', fontsize)
xticks(0:20:100)
xtickangle(0)

% load(fullfile(base, 'model_comparisons', 'Pop_Rate_Timbre_Best.mat'), ...
% 	"pop_rate_timbre_best")
% nneurons = [1:4 5:5:50];
% accuracy_bad = [pop_rate_timbre_best(15:28).accuracy];
% plot(nneurons, accuracy_bad);
% 
% hold on 
% accuracy_good = [pop_rate_timbre_best(1:14).accuracy];
% plot(nneurons, accuracy_good);
% xlabel('# Neurons in Model')
% ylabel('Model Accuracy')
% grid on
% box off
% hleg = legend('Worst', 'Best', 'fontsize', legsize, 'location', 'best', 'box', 'off');
% hleg.ItemTokenSize = [8, 8];
% set(gca, 'fontsize', fontsize)



%% I. Subset MTF groups 
% h(10) = subplot(4, 3, 11);
% 
% load(fullfile(base, 'model_comparisons', 'Pop_Rate_Timbre_Subset_MTF.mat'), ...
% 	"accuracy_all", "MTF_names")
% 
% hold on
% for iCF = 1:5
% 	if iCF == 1
% 		RGB = hex2rgb(colorsTimbre);
% 	else
% 		RGB = hex2rgb(colorsMTF{iCF-1});
% 	end
% 	swarmchart(ones(size(accuracy_all, 2), 1)*iCF, accuracy_all(iCF,:), scattersize-4, RGB)
% end
% xlabel('MTF Groups')
% xticks(1:5)
% xticklabels(MTF_names)
% ylim([0.78 1])
% max_acc = max(accuracy_all, [], 2);
% plot(1:5, max_acc, 'k')
% mean_acc = median(accuracy_all, 2);
% plot(1:5, mean_acc, 'k')
% set(gca, 'fontsize', fontsize)
% ylabel('Model Accuracy')
% grid on



%% Set position

left = linspace(0.1, 0.74, 3);
bottom = [0.07 0.36 0.64];
height1 = 0.21;
height = 0.32;
width = 0.23;

set(h(3), 'Position', [left(1) bottom(3) 0.87 height/3])
set(h(2), 'Position', [left(1) bottom(3)+height/3 0.87 height/3])
set(h(1), 'Position', [left(1) bottom(3)+height/3*2 0.87 height/3])

set(h(4), 'position', [left(1) bottom(2) width height1])
set(h(5), 'position', [left(1) bottom(1) width height1])

set(h(6), 'position', [left(2) bottom(2) width height1])
set(h(7), 'position', [left(3) bottom(2) width height1])
%set(h(8), 'position', [left(3) bottom(2) width height])
set(h(9), 'position', [left(2) bottom(1) width height1])
set(h(10), 'position', [left(3) bottom(1) width height1])
%set(h(11), 'position', [left(3) bottom(1) width height])


%% Annotate 

left = linspace(0.03, 0.68, 3);
bottom = [0.3, 0.57, 0.97];

label = {'A', 'B', 'C'};
for ii = 1
	annotation('textbox',[left(ii) bottom(3) 0.0826 0.0385],'String',label{ii},...
		'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
end
label = {'B', 'C', 'D'};
for ii = 1:3
	annotation('textbox',[left(ii) bottom(2) 0.0826 0.0385],'String',label{ii},...
		'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
end
label = {'E', 'F', 'G'};
for ii = 1:3
	annotation('textbox',[left(ii) bottom(1) 0.0826 0.0385],'String',label{ii},...
		'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
end

%% Save figure 

if save_fig == 1
	filename = 'fig5_timbre_population_rate';
	save_figure(filename)
end
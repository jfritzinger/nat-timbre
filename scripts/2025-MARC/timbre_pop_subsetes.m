%% timbre_pop_subsets
clear
save_fig = 1;

%% Load in data

[base, ~, ~, ppi] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Pop_Rate_Timbre_Subset_MTF2.mat'), ...
	"accuracy", "std_acc", "MTF_names")
acc_MTF = accuracy;
std_MTF = std_acc;
load(fullfile(base, 'model_comparisons', 'Pop_Rate_Timbre_Subset_CF2.mat'), ...
	"accuracy", "std_acc", "CF_groups")


%% Set up figure

figure('Position',[50, 50, 10.5*ppi, 3.5*ppi])
tiledlayout(1, 2, 'Padding','compact')
linewidth = 2;
fontsize = 18;
titlesize = 20;
legsize = 20;
scattersize = 40;
colorsData = {"#0072BD", "#D95319"};
colorsTimbre = '#1b9e77';
colorsCF = [27/256, 158/256, 119/256;0.0660 0.4430 0.7450; 0.8660 0.3290 0.0000;0.9290 0.6940 0.1250];

%% Plot CFs 

num_subset = [1, 2, 3, 4, 5:5:90];
nexttile
hold on
for i = 1:4
	plot(num_subset, accuracy(i,:), 'linewidth', 2, 'Color',colorsCF(i,:))
end
ylabel('Accuracy')
xlabel('# Neurons in Model')
legend('Best Neurons', 'CF = 0-2 kHz', 'CF = 2-4 kHz', 'CF = 4-8 kHz', 'Location','southeast')
grid on 
box off
set(gca, 'fontsize', fontsize)
xticks(0:20:100)
xtickangle(0)

%% Plot MTFs  

colorsMTF = {'#1b9e77', '#648FFF', '#DC267F', '#785EF0', '#FFB000'};

nexttile
hold on
for iMTF = 1:5
	plot(num_subset, acc_MTF(iMTF,:), 'Color', colorsMTF{iMTF}, 'linewidth', 2)
end
ylabel('Accuracy')
xlabel('# Neurons in Model')
legend({'Best Neurons', MTF_names{2:5}}, 'Location','southeast')
grid on 
box off
set(gca, 'fontsize', fontsize)
xticks(0:20:100)
xtickangle(0)


%% Annotate / position



%% Save figure 

if save_fig == 1
	filename = 'timbre_pop_subsets';
	save_figure_MARC(filename)
end

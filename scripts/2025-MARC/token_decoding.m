%% token_decoding
clear
save_fig = 1;

%% Load in data 

[base, ~, ~, ppi] = getPathsNT;
load(fullfile(base, 'model_comparisons', 'Pop_Rate.mat'), "pop_rate")
load(fullfile(base, 'model_comparisons', 'Pop_Time.mat'), ...
	"accur_all","C_all", "num_neurons", "mean_acc", "std_acc")

%% Set up figure

figure('Position',[50, 50, 10.5*ppi, 4*ppi])
%tiledlayout(1, 2, 'Padding','compact')
linewidth = 2;
fontsize = 18;
titlesize = 20;
legsize = 20;
scattersize = 40;
colorsData = {"#0072BD", "#D95319"};
colorsTimbre = '#1b9e77';

% Calculate accuracy
accuracy(1) = sum(diag(pop_rate.C)) / sum(pop_rate.C(:)); % Calculate accuracy
x = 1:75;
h(1) = subplot(1, 3, 1);
pcolor(x, x, pop_rate.C, 'EdgeColor','none')
clim([0 10])
a=colorbar;
a.Label.String = '# Accurate Predictions';
set(gca, 'FontSize', fontsize)
% colormap(brewermap([],"GnBu"))
% colormap(brewermap([],"BuGn"))
colormap(h(1), brewermap([],"Blues"))
title('Rate Model')

% Annotations
hold on
x = [1, 41, 41, 1, 1];
y = [1, 1, 41, 41, 1];
plot(x, y, 'k-', 'LineWidth', 1);

x = [41, 75, 75, 41, 41];
y = [41, 41, 75, 75, 41];
plot(x, y, 'k-', 'LineWidth', 1);
xticklabels([])
yticklabels([])

% Plot timing results 
h(2) = subplot(1, 3, 2);
plot(num_neurons, accur_all, 'LineWidth',2)
ylabel('Overall Accuracy')
xlabel('                               # Neurons in Model')
title('                               Timing Model')
set(gca, 'fontsize', fontsize)
grid on
ylim([ 0 1])

h(3) = subplot(1, 3, 3);
for ind = 1:10
	C_diag(ind,:) = diag(C_all{ind});
end
pcolor(num_neurons, 1:75,C_diag', 'EdgeColor','none')
hold on
yline(41, 'k', 'LineWidth',1)
colormap(h(3), "parula")
colorbar
set(gca, 'fontsize', fontsize)
yticklabels([])

% Arrange 
left = [0.08 0.55 0.78];
bottom = 0.16;
height = 0.75;
width = 0.15;

set(h(1), 'position', [left(1) bottom 0.28 height])
set(h(2), 'position', [left(2) bottom width height])
set(h(3), 'position', [left(3) bottom width height])

annotation('textbox',[0.1 0.03 0.25 0.058],...
	'String','Bassoon','FontSize',18,...
	'EdgeColor','none', 'FontWeight','bold');
annotation('textbox',[0.27 0.03 0.25 0.058],...
	'String','Oboe','FontSize',18,...
	'EdgeColor','none', 'FontWeight','bold');
annotation('textbox',[0.03 0.7 0.25 0.058],...
	'String','Oboe','FontSize',18,...
	'EdgeColor','none', 'Rotation',90, 'FontWeight','bold');
annotation('textbox',[0.03 0.21 0.25 0.058],...
	'String','Bassoon','FontSize',18,...
	'EdgeColor','none', 'Rotation',90, 'FontWeight','bold');
% annotation("textarrow", [0.11 0.3532], [0.122 0.122],...
% 	"String", "F0",'FontSize',18)

annotation("textarrow", [0.26 0.3532], [0.122 0.122],...
	"String", "F0",'FontSize',18)
annotation("textarrow", [0.11 0.225], [0.122 0.122],...
	"String", "F0",'FontSize',18)
annotation("textarrow", [0.06 0.06], [0.22 0.55],...
	"String", "F0",'FontSize',18)
annotation("textarrow", [0.06 0.06], [0.63 0.9],...
	"String", "F0",'FontSize',18)


annotation('textbox',[0.74 0.66 0.25 0.058],...
	'String','Oboe','FontSize',18,...
	'EdgeColor','none', 'Rotation',90, 'FontWeight','bold');
annotation('textbox',[0.74 0.21 0.25 0.058],...
	'String','Bassoon','FontSize',18,...
	'EdgeColor','none', 'Rotation',90, 'FontWeight','bold');
annotation("textarrow", [0.76 0.76], [0.22 0.55],...
	"String", "F0",'FontSize',18)
annotation("textarrow", [0.76 0.76], [0.63 0.9],...
	"String", "F0",'FontSize',18)


%% Save figure

if save_fig == 1
	filename = 'token_decoding';
	save_figure_MARC(filename)
end



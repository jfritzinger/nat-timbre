%% token_decoding
clear
save_fig = 1;

%% Load in data 

[base, ~, ~, ppi] = getPathsNT;
load(fullfile(base, 'model_comparisons', 'Pop_Rate.mat'), "pop_rate")
load(fullfile(base, 'model_comparisons', 'Pop_Time.mat'), ...
	"accur_all","C_all", "num_neurons", "mean_acc", "std_acc")

%% Figure

figure('position', [50, 50, 4.7*ppi, 4*ppi])
linewidth = 2;
fontsize = 18;
titlesize = 20;
legsize = 20;
scattersize = 40;

% Seperate: 

x = 1:75;
pcolor(x, x, C_all{end}, 'EdgeColor','none')
clim([0 3])
a=colorbar;
a.Label.String = '# Accurate Predictions';
set(gca, 'FontSize', fontsize)
% colormap(brewermap([],"GnBu"))
% colormap(brewermap([],"BuGn"))
colormap(brewermap([],"Blues"))
title('Timing Model')


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

annotation("ellipse", [0.3741 0.5712 0.07805 0.3341])
annotation("ellipse", [0.4917 0.414 0.2854 0.1227])
%%
annotation('textbox',[0.17 0.03 0.25 0.058],...
	'String','Bassoon','FontSize',18,...
	'EdgeColor','none', 'FontWeight','bold');
annotation('textbox',[0.56 0.03 0.25 0.058],...
	'String','Oboe','FontSize',18,...
	'EdgeColor','none', 'FontWeight','bold');
annotation('textbox',[0.03 0.65 0.25 0.058],...
	'String','Oboe','FontSize',18,...
	'EdgeColor','none', 'Rotation',90, 'FontWeight','bold');
annotation('textbox',[0.03 0.18 0.25 0.058],...
	'String','Bassoon','FontSize',18,...
	'EdgeColor','none', 'Rotation',90, 'FontWeight','bold');


annotation("textarrow", [0.55 0.77], [0.08 0.08],...
	"String", "F0",'FontSize',18)
annotation("textarrow", [0.17 0.47], [0.08 0.08],...
	"String", "F0",'FontSize',18)
annotation("textarrow", [0.08 0.08], [0.18 0.52],...
	"String", "F0",'FontSize',18)
annotation("textarrow", [0.08 0.08], [0.6 0.9],...
	"String", "F0",'FontSize',18)

%% Save figure

if save_fig == 1
	filename = 'token_decoding2';
	save_figure_MARC(filename)
end


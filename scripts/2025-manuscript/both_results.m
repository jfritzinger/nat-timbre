%% both_results
clear
save_fig = 1;

%% Load in data

[base, ~, ~, ppi] = getPathsNT;
load(fullfile(base, 'model_comparisons', 'Pop_Rate2.mat'), "pop_rate")
load(fullfile(base, 'model_comparisons', 'Pop_Time.mat'), ...
	"accur_all","C_all", "num_neurons", "mean_acc", "std_acc")

load(fullfile(base, 'model_comparisons', ['Neuron_Time_F0_' ...
	'Oboe' '.mat']), "neuron_time_F0")
neuron_time_F0_oboe = neuron_time_F0;
load(fullfile(base, 'model_comparisons', ['Neuron_Rate_F0_' ...
	'Oboe' '.mat']), "neuron_rate_F0")
neuron_rate_F0_oboe = neuron_rate_F0;

load(fullfile(base, 'model_comparisons', ['Neuron_Rate_F0_' ...
	'Bassoon' '.mat']), "neuron_rate_F0")
load(fullfile(base, 'model_comparisons', ['Neuron_Time_F0_' ...
	'Bassoon' '.mat']), "neuron_time_F0")

load(fullfile(base, 'model_comparisons', 'Neuron_Rate_Timbre_All.mat'), ...
	"neuron_rate_timbre")
load(fullfile(base, 'model_comparisons', 'Neuron_Time_Timbre_All.mat'), ...
	"neuron_time_timbre")

%% Set up figure

figure('Position',[50, 50, ppi*6.5, ppi*2.5])
scattersize = 5;
titlesize = 9;
fontsize = 8;
labelsize = 12;
legsize = 7;
colorsData = {"#0072BD", "#D95319"};
colorsMTF = {'#648FFF', '#DC267F', '#785EF0', '#FFB000'};
colorsTimbre = '#1b9e77';
linewidth = 1;

% Single-unit analysis

target = {'Bassoon', 'Oboe'};
h_ind = [0, 2];
h_ind2 = [1, 2; 8, 9];
for itype = 1:2
	for itarget = 1:2

		% Load data
		if itype == 1
			if itarget == 1
				neuron_timbre = neuron_rate_timbre;
				neuron_F0 = neuron_rate_F0;
			else
				neuron_timbre = neuron_rate_timbre;
				neuron_F0 = neuron_rate_F0_oboe;
			end
		else
			if itarget == 1
				neuron_timbre = neuron_time_timbre;
				neuron_F0 = neuron_time_F0;
			else
				neuron_timbre = neuron_time_timbre;
				neuron_F0 = neuron_time_F0_oboe;
			end
		end

		% Simplify neuron_rate_F0 to match timbre
		putatives_timbre = {neuron_timbre.putative};
		putatives_F0 = {neuron_F0.putative};
		for i = 1:numel(putatives_timbre)
			idx = find(strcmp(putatives_F0, putatives_timbre{i}));
			if ~isempty(idx)
				indices(i) = idx;
			else
				indices(i) = 0; % or NaN, to indicate not found
			end
		end
		indices(indices==0) = [];
		accuracy_timbre = [neuron_timbre.accuracy];
		accuracy_F0 = [neuron_F0(indices).accuracy];

		h(h_ind(itype)+itarget) = subplot(2, 7, h_ind2(itype,itarget));
		scatter(accuracy_timbre, accuracy_F0, scattersize, 'filled', ...
			'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.5);
		if itarget == 1 && itype == 2
			xlabel('                   Instrument Accuracy')
		end
		if itarget == 1 && itype == 1
			ylabel('F0 Accuracy                                ')
		end
		hold on
		xlim([0.4 1])
		ylim([0 1])
		yticks(0:0.2:1)
		xticks(0:0.2:1)
		if itype == 1
			xticklabels([])
		end
		if itarget == 2
			yticklabels([])
		end
		xtickangle(0)
		grid on

		if itype == 1
			title(target{itarget})
		end

		mdl = fitlm(accuracy_timbre, accuracy_F0);
		x = linspace(0, 1, 20);
		y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
		plot(x, y, ':k', 'linewidth', 1)
		hleg = legend('', ...
			sprintf('p=%0.04f',mdl.Coefficients{2,4}), 'fontsize', legsize, ...
			'location', 'northwest', 'box', 'off');
		hleg.ItemTokenSize = [12, 8];
		set(gca, 'fontsize', fontsize)

	end
end

% Rate Population

% Calculate accuracy
accuracy(1) = sum(diag(pop_rate.C)) / sum(pop_rate.C(:)); % Calculate accuracy
x = 1:75;
h(5) = subplot(2, 7, [3, 4, 10, 11]);
pcolor(x, x, pop_rate.C, 'EdgeColor','none')
clim([0 3])
%a=colorbar;
%a.Label.String = '# Accurate Predictions';
set(gca, 'FontSize', fontsize)
% colormap(brewermap([],"GnBu"))
% colormap(brewermap([],"BuGn"))
colormap(h(1), brewermap([],"Blues"))
title('Rate Model')
xlabel('Actual F0 (Hz)')
ylabel('Predicted F0 (Hz)')

% Annotations
hold on
x = [1, 41, 41, 1, 1];
y = [1, 1, 41, 41, 1];
plot(x, y, 'k-', 'LineWidth', 1);
x = [41, 75, 75, 41, 41];
y = [41, 41, 75, 75, 41];
plot(x, y, 'k-', 'LineWidth', 1);
%xticklabels([])
%yticklabels([])
F0s_b = getF0s('Bassoon');
F0s_o = getF0s('Oboe');
F0s_all = round([F0s_b; F0s_o]);
xticks([1 20 35 41 60 75])
xticklabels([F0s_all(1) F0s_all(20) F0s_all(35) F0s_all(41) F0s_all(60) F0s_all(75)])
yticks([1 20 35 41 60 75])
yticklabels([F0s_all(1) F0s_all(20) F0s_all(35) F0s_all(41) F0s_all(60) F0s_all(75)])

% Timing Population

% Plot timing results
h(6) = subplot(2, 7, [5, 6, 12, 13]);
x = 1:75;
pcolor(x, x, C_all{end}, 'EdgeColor','none')
accuracy(2) = sum(diag(C_all{end})) / sum(C_all{end}, "all"); % Calculate accuracy

clim([0 3])
a=colorbar;
a.Label.String = '# Accurate Predictions';
set(gca, 'FontSize', fontsize)
% colormap(brewermap([],"GnBu"))
% colormap(brewermap([],"BuGn"))
colormap(brewermap([],"Blues"))
title('Timing Model')
xlabel('Actual F0 (Hz)')


hold on
x = [1, 41, 41, 1, 1];
y = [1, 1, 41, 41, 1];
plot(x, y, 'k-', 'LineWidth', 1);

x = [41, 75, 75, 41, 41];
y = [41, 41, 75, 75, 41];
plot(x, y, 'k-', 'LineWidth', 1);
yticklabels([])
xticks([1 20 35 45 60 75])
xticklabels([F0s_all(1) F0s_all(20) F0s_all(35) F0s_all(45) F0s_all(60) F0s_all(75)])

% % Timing F0s
% h(7) = subplot(2, 7, [7, 14]);
% for ind = 1:10
% 	C_diag(ind,:) = diag(C_all{ind});
% end
% pcolor(num_neurons, 1:75,C_diag', 'EdgeColor','none')
% hold on
% yline(41, 'k', 'LineWidth',1)
% colormap(h(7), "parula")
% c = colorbar;
% c.Label.String = '# Accurate Predictions';
% set(gca, 'fontsize', fontsize)
% yticklabels([])
% xlabel('# Neurons')

% Arrange figure

left = [0.09 0.23 0.46 0.69];
bottom = [0.22 0.58];
height = 0.33;
height2 = 0.68;
width = 0.12;
width2 = 0.21;

set(h(1), 'position', [left(1) bottom(2) width height])
set(h(2), 'position', [left(2) bottom(2) width height])
set(h(3), 'position', [left(1) bottom(1) width height])
set(h(4), 'position', [left(2) bottom(1) width height])

set(h(5), 'position', [left(3) bottom(1) width2 height2])
set(h(6), 'position', [left(4) bottom(1) width2 height2])
% set(h(7), 'position', [left(5) bottom(1) width-0.04 height2])
%
annotation('textbox',[left(3)+0.01 0.03 0.25 0.058],...
	'String','Bassoon','FontSize',titlesize,...
	'EdgeColor','none', 'FontWeight','bold');
annotation('textbox',[left(3)+0.14 0.03 0.25 0.058],...
	'String','Oboe','FontSize',titlesize,...
	'EdgeColor','none', 'FontWeight','bold');
annotation('textbox',[left(4)+0.01 0.03 0.25 0.058],...
	'String','Bassoon','FontSize',titlesize,...
	'EdgeColor','none', 'FontWeight','bold');
annotation('textbox',[left(4)+0.14 0.03 0.25 0.058],...
	'String','Oboe','FontSize',titlesize,...
	'EdgeColor','none', 'FontWeight','bold');

annotation('textbox',[left(3)-0.075 0.64 0.25 0.058],...
	'String','Oboe','FontSize',titlesize,...
	'EdgeColor','none', 'Rotation',90, 'FontWeight','bold');
annotation('textbox',[left(3)-0.075 0.21 0.25 0.058],...
	'String','Bassoon','FontSize',titlesize,...
	'EdgeColor','none', 'Rotation',90, 'FontWeight','bold');
% annotation('textbox',[left(5)-0.015 0.58 0.25 0.058],...
% 	'String','Oboe F0','FontSize',fontsize,...
% 	'EdgeColor','none', 'Rotation',90, 'FontWeight','bold');
% annotation('textbox',[left(5)-0.015 0.18 0.25 0.058],...
% 	'String','Bassoon F0','FontSize',fontsize,...
% 	'EdgeColor','none', 'Rotation',90, 'FontWeight','bold');

% annotation("textarrow", [left(3)+0.02 left(3)+0.1], [0.122 0.122],...
% 	"String", "F0",'FontSize',fontsize, 'HeadLength',7, 'HeadWidth',7)
% annotation("textarrow", [left(3)+0.13 left(3)+0.2], [0.122 0.122],...
% 	"String", "F0",'FontSize',fontsize, 'HeadLength',7, 'HeadWidth',7)
% annotation("textarrow", [left(4)+0.02 left(4)+0.1], [0.122 0.122],...
% 	"String", "F0",'FontSize',fontsize, 'HeadLength',7, 'HeadWidth',7)
% annotation("textarrow", [left(4)+0.13 left(4)+0.2], [0.122 0.122],...
% 	"String", "F0",'FontSize',fontsize, 'HeadLength',7, 'HeadWidth',7)
% annotation("textarrow", [left(3)-0.015  left(3)-0.015], [0.22 0.55],...
% 	"String", "F0",'FontSize',fontsize, 'HeadLength',7, 'HeadWidth',7)
% annotation("textarrow", [left(3)-0.015 left(3)-0.015], [0.63 0.9],...
% 	"String", "F0",'FontSize',fontsize, 'HeadLength',7, 'HeadWidth',7)

annotation('textbox',[left(1)-0.06 0.68 0.25 0.058],...
	'String','Rate','FontSize',titlesize, 'Rotation',90,...
	'EdgeColor','none', 'FontWeight','bold');
annotation('textbox',[left(1)-0.06 0.24 0.25 0.058],...
	'String','Timing','FontSize',titlesize,'Rotation',90,...
	'EdgeColor','none', 'FontWeight','bold');

% Annotate 
annotation('textbox',[left(1)-0.08 0.98 0.0826 0.0385],'String','A',...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[left(3)-0.08 0.98 0.0826 0.0385],'String','B',...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[left(4)-0.03 0.98 0.0826 0.0385],'String','C',...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
% annotation('textbox',[left(5)-0.04 0.98 0.0826 0.0385],'String','D',...
% 	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');

%% Save figure

if save_fig == 1
	filename = 'fig10_both_results';
	save_figure(filename)
end
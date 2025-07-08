%%pitch_single_unit 
% clear
% save_fig = 1;
% 
% %% Load in data 
% 
% [base, datapath, ~, ppi] = getPathsNT();
% 
% % Load in oboe
% load(fullfile(base, 'model_comparisons', 'Neuron_Time_F0_Oboe.mat'), ...
% 	"neuron_time_F0")
% neuron_time_F0_oboe = neuron_time_F0;
% load(fullfile(base, 'model_comparisons','Neuron_Rate_F0_Oboe.mat'), ...
% 	"neuron_rate_F0")
% neuron_rate_F0_oboe = neuron_rate_F0;
% 
% % Load in bassoon
% load(fullfile(base, 'model_comparisons', 'Neuron_Time_F0_Bassoon.mat'), ...
% 	"neuron_time_F0")
% load(fullfile(base, 'model_comparisons','Neuron_Rate_F0_Bassoon.mat'), ...
% 	"neuron_rate_F0")
% load(fullfile(base,'model_comparisons', 'Data_NT_3.mat'), 'nat_data')
% 
% accuracy_rate = [neuron_rate_F0.accuracy]*100;
% accuracy_time = [neuron_time_F0.accuracy]*100;
% accuracy_rate_oboe = [neuron_rate_F0_oboe.accuracy]*100;
% accuracy_time_oboe = [neuron_time_F0_oboe.accuracy]*100;
% CFs = [neuron_time_F0.CF];
% MTFs = {neuron_time_F0.MTF};
% F0s = getF0s('Bassoon');
% 
% % Analysis 
% sessions = readtable(fullfile(base, 'Data_Table.xlsx'), ...
% 	'PreserveVariableNames',true);
% [ac, ind_high] = sort(accuracy_time, 'ascend');
% for ii = 1:length(accuracy_time)
% 
% 	putative = neuron_time_F0(ind_high(ii)).putative;
% 	load(fullfile(datapath, [putative '.mat']), 'data');
% 	s_ind = strcmp(sessions.Putative_Units, putative);
% 	CF = sessions.CF(s_ind);
% 	params_NT = data(7, 2);
% 	data_NT = analyzeNT(params_NT{1});
% 	temporal = analyzeNT_Temporal(data_NT, CF);
% 	r_splithalf(:, ii) = temporal.r_splithalf;
% 	VS_all2(:, ii) = temporal.VS;
% end

%% Create figure

figure('Position',[50,600,8*ppi,4.75*ppi])
linewidth = 1;
fontsize = 8;
legsize = 7;
labelsize = 12; 
colorsData = {"#0072BD", "#D95319"};
colorsMTF = {'#648FFF', '#DC267F', '#785EF0', '#FFB000'};
colorsPitch = '#0072B2';
scattersize = 10;

%% A. Bassoon rate vs timing 

edges = linspace(0, 100, 101);
h(1) = subplot(4, 4, 1);
set(gca, 'fontsize', fontsize)
histogram(accuracy_rate,edges, 'FaceColor', colorsPitch)
xlabel('Rate Accuracy (%)')
grid on
box off 
set(gca, 'fontsize', fontsize)
xlim([0 25])

h(2) = subplot(4, 4, 2);
set(gca, 'fontsize', fontsize)
histogram(accuracy_time, edges, 'Orientation','horizontal', 'FaceColor', colorsPitch)
ylabel('Timing Accuracy (%)')
grid on
box off
set(gca, 'fontsize', fontsize)
ylim([0 60])

h(3) = subplot(4, 4, 3);
hold on
scatter(accuracy_rate, accuracy_time, scattersize, 'filled',...
	'MarkerEdgeColor','k', 'MarkerFaceColor', colorsPitch, 'MarkerFaceAlpha',0.5)
plot([0 100], [0 100], 'k')
grid on
xlim([0 25])
ylim([0 60])
set(gca, 'fontsize', fontsize)
xticklabels([])
yticklabels([])
title('Bassoon')

fprintf('Bassoon, Rate: Mean = %0.02f, Max = %0.02f\n', ...
	mean(accuracy_rate), max(accuracy_rate))
fprintf('Bassoon, Timing: Mean = %0.02f, Max = %0.02f\n', ...
	mean(accuracy_time), max(accuracy_time))

%% B. Oboe rate vs timing 

h(4) = subplot(4, 4, 4);
edges = linspace(0, 100, 101);
chance = 1/length(neuron_rate_F0_oboe(1).rate)*100;
hold on
accuracy = accuracy_rate_oboe;
histogram(accuracy,edges)
xline(chance, 'k', 'LineWidth',linewidth)
xline(mean(accuracy), 'r', 'LineWidth',linewidth)
xlim([0 25])
xlabel('Rate Accuracy (%)')
grid on
set(gca, 'fontsize', fontsize)
fprintf('Oboe, Rate: Mean = %0.02f, Max = %0.02f\n', ...
	mean(accuracy), max(accuracy))

h(5) = subplot(4, 4, 5);
edges = linspace(0, 100, 101);
chance = 1/length(neuron_rate_F0_oboe(1).rate)*100;
hold on
accuracy = accuracy_time_oboe;
histogram(accuracy,edges, 'Orientation','horizontal')
yline(chance, 'k', 'LineWidth',linewidth)
yline(mean(accuracy), 'r', 'LineWidth',linewidth)
ylim([0 60])
ylabel('Timing Accuracy (%)')
grid on
set(gca, 'fontsize', fontsize)


h(6) = subplot(4, 4, 6);
scatter(accuracy_rate_oboe, accuracy_time_oboe, scattersize, 'filled', 'MarkerEdgeColor','k', ...
	'MarkerFaceAlpha',0.5)
hold on
plot([0, 100], [0, 100], 'k')
xlim([0 25])
ylim([0 60])

yticks(0:5:100)
xticks(0:5:100)
title('Oboe')
set(gca, 'fontsize', fontsize)
xticklabels([])
yticklabels([])
grid on


fprintf('Oboe, Timing: Mean = %0.02f, Max = %0.02f\n', ...
	mean(accuracy), max(accuracy))

%% C. CFs

h(7) = subplot(4, 4, 7);

scatter(CFs, accuracy_time, scattersize, 'filled', 'MarkerEdgeColor',...
	'k', 'MarkerFaceColor', colorsPitch, 'MarkerFaceAlpha',0.5);
set(gca, 'xscale', 'log')
hold on
set(gca, 'fontsize', fontsize)
xticks([200 500 1000 2000 5000 10000])
xticklabels([0.2 0.5 1 2 5 10])
xlabel('CF (kHz)')
ylabel('Bassoon Accuracy (%)')
grid on

%% D. MTFs

h(8) = subplot(4, 4, 8);

MTF_types = unique(MTFs);
hold on
for iMTF = 1:4
	if iMTF == 4
		ind = strcmp(MTFs, MTF_types{iMTF}) | strcmp(MTFs, MTF_types{iMTF+1});
	else
		ind = strcmp(MTFs, MTF_types{iMTF});
	end
	accur = accuracy_time(ind);
	num_units = length(accur);

	RGB = hex2rgb(colorsMTF{iMTF});
	swarmchart(ones(num_units, 1)*iMTF, accur, scattersize, RGB)

	mean_vals(iMTF) = mean(accur);
	std_vals(iMTF) = std(accur)/sqrt(length(accur));
end
xticks(1:4)
xticklabels({'BE', 'BS', 'F', 'H'})
xlabel('MTF Groups')
ylabel('Bassoon Accuracy (%)')

% tableMTF = table(MTFs', accuracy(1,:)');
% anova(tableMTF, 'Var2')
% [~,~,stats] = anova1(accuracy(1,:), MTFs);
% [c,~,~,gnames] = multcompare(stats);
grid on
set(gca, 'fontsize', fontsize)

%% E.  RVF 

h(9) = subplot(4, 4, 9);
VS_all = zeros(3,1);
for ii = 1:297
	VS_1 = nat_data(ii).bass_VS;
	if ~isempty(VS_1)
		PC2_score(ii,:) = nat_data(ii).RVF_PC2;
	else
		PC2_score(ii,:) = 0;
	end
end
PC2_score(PC2_score==0) = [];

% Plot PCA2 score vs beta weights 
scatter(PC2_score, accuracy_time, scattersize, 'filled',...
	'MarkerEdgeColor','k', 'MarkerFaceColor', colorsPitch)
xlabel('PC2 RVF score')
ylabel('Bassoon Accuracy (%)')
hold on

mdl = fitlm(PC2_score, accuracy_time);
x = linspace(-2, 1, 20);
y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
plot(x, y, ':k')
% hleg = legend('Neuron', ...
% 	sprintf('p=%0.04f',mdl.Coefficients{2,4}), 'fontsize', 7, ...
% 	'location', 'northwest', 'box', 'off');
% hleg.ItemTokenSize = [8, 8];
grid on
set(gca, 'fontsize', fontsize)
grid on



%% F. Imagesc accuracy

h(10) = subplot(4, 4, 10);

[~,best_ind] = sort(abs(accuracy_time), 'ascend' );
for ii = 1:287
	C_acc(ii,:) = diag(neuron_time_F0(best_ind(ii)).C);
end

x = getF0s('Bassoon');
y = 1:287;

pcolor(x, y, C_acc, 'EdgeColor','none', 'EdgeAlpha',0)
set(gca, 'xscale', 'log')
xticks([60 100 200 350 550])
xlabel('F0 (Hz)')
ylabel('Neurons Sorted by Accuracy')
%c = colorbar;
title('Bassoon Accuracy',...
	'HorizontalAlignment','center')
set(gca, 'fontsize', fontsize)
%c.Label.String = '# Accurate Predictions';

% Analysis
best = mean(C_acc(end-15:end,:));


%% Comparison of oboe and bassoon

h(11) = subplot(4, 4, 11);
[sesh, num_data] = getF0Sessions(nat_data, 'Invariant');
putative = {nat_data(sesh).putative};
nneurons = length(sesh);

for i = 1:nneurons

	putative1 = putative{i};
	ind_bass = find(strcmp(putative1, {neuron_time_F0.putative}));
	ind_oboe = find(strcmp(putative1, {neuron_time_F0_oboe.putative}));

	accuracy_oboe(i) = neuron_time_F0_oboe(ind_oboe).accuracy*100;
	accuracy_bass1(i) = neuron_time_F0(ind_bass).accuracy*100;
end

scatter(accuracy_bass1, accuracy_oboe, scattersize, 'filled', 'MarkerEdgeColor','k', ...
	'MarkerFaceAlpha',0.5)
hold on
plot([0, 100], [0, 100], 'k')
xlim([0 65])
ylim([0 65])
set(gca, 'fontsize', fontsize)
grid on
ylabel('Oboe Accuracy (%)')
xlabel('Bassoon Accuracy (%)')

%% Imagesc reliability 

h(12) = subplot(4, 4, 12);
pcolor(F0s, 1:length(accuracy_time), r_splithalf', 'EdgeColor','none', 'EdgeAlpha',0)
set(gca, 'xscale', 'log')
xticks([60 100 200 350 550])
xlabel('F0 (Hz)')
title('Reliability')
yticklabels([])
set(gca, 'fontsize', fontsize)

%% Reliability 

h(13) = subplot(4, 4, 13);
mean_splithalf = mean(r_splithalf, 1);
scatter(ac, mean_splithalf, scattersize, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.5)
ylabel('Average Reliability')
xlabel('Accuracy')
set(gca, 'fontsize', fontsize)
grid on
ylim([0 0.6])

%% Imagesc vector strength 

h(14) = subplot(4, 4, 14);
pcolor(F0s, 1:length(accuracy_time), VS_all2', 'EdgeColor','none', 'EdgeAlpha',0)
c = colorbar;
c.Label.String = 'Accuracy / Reliability / VS';
set(gca, 'xscale', 'log')
xticks([60 100 200 350 550])
xlabel('F0 (Hz)')
title('Vector Strength')
set(gca, 'fontsize', fontsize)
grid on
yticklabels([])


%% Vector strength

h(15) = subplot(4, 4, 15);
mean_VS = mean(VS_all2, 1);
scatter(ac, mean_VS, scattersize, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.5)
hold on
% mdl = fitlm(mean_VS, accuracy_time);
% x = linspace(0, 0.6, 10);
% y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
% plot(x, y, ':k')

ylabel('Average VS')
xlabel('Accuracy')
set(gca, 'fontsize', fontsize)
grid on
ylim([0 0.6])

%% Arrange 

left = [0.095 0.325 0.52 0.66 0.8];
bottom = linspace(0.08, 0.73, 3);
height = 0.22;
width = 0.12;

fig_position = [left(1),0.65,0.16,0.3];
nb_position = [fig_position(1),fig_position(2)-0.08,fig_position(3),0.06];
wb_position = [fig_position(1)-0.055,fig_position(2),0.045,fig_position(4)];
set(h(3), 'Position', fig_position)
set(h(1), 'Position', nb_position)
set(h(2), 'Position', wb_position)

fig_position = [left(1),bottom(1)+0.08,0.16,0.3];
nb_position = [fig_position(1),fig_position(2)-0.08,fig_position(3),0.06];
wb_position = [fig_position(1)-0.055,fig_position(2),0.045,fig_position(4)];
set(h(6), 'Position', fig_position)
set(h(4), 'Position', nb_position)
set(h(5), 'Position', wb_position)


set(h(7), 'position', [left(2) bottom(3) width height])
set(h(8), 'position', [left(2) bottom(2) width height])
set(h(9), 'position', [left(2) bottom(1) width height])
set(h(10), 'position', [left(3) bottom(2) width height*2+0.1]) % image
set(h(11), 'position', [left(3) bottom(1) width-0.02 height])
set(h(12), 'position', [left(4) bottom(2) width height*2+0.1]) % image
set(h(13), 'position', [left(4)+0.02 bottom(1) width-0.02 height])
set(h(14), 'position', [left(5) bottom(2) width height*2+0.1]) % image
set(h(15), 'position', [left(5)+0.04 bottom(1) width-0.02 height])

%% Set labels 

labelleft= [0 0.26 0.48 0.64 0.8];
labelbottom = linspace(0.29, 0.94, 3);
annotation('textbox',[labelleft(1) labelbottom(3) 0.071 0.058],...
	'String','A','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(1) 0.47 0.071 0.058],...
	'String','B','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(2) labelbottom(3) 0.071 0.058],...
	'String','C','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(2) labelbottom(2) 0.071 0.058],...
	'String','D','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(2) labelbottom(1) 0.071 0.058],...
	'String','E','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(3) labelbottom(3) 0.071 0.058],...
	'String','F','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');

annotation('textbox',[labelleft(3) labelbottom(1) 0.071 0.058],...
	'String','G','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(4) labelbottom(1) 0.071 0.058],...
	'String','H','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(5) labelbottom(1) 0.071 0.058],...
	'String','I','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');


%% Save figure 

if save_fig == 1
	filename = 'fig6_single_neuron_F0';
	save_figure(filename)
end
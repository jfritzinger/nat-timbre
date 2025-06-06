%% pitch_single_unit_oboe
clear
save_fig = 1;

%% Get paths

[base, ~, savepath, ppi] = getPathsNT();

%% Load in data 

load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')

filename = 'Neuron_Time_F0_Bassoon';
filepath = fullfile(base, 'model_comparisons',filename);
load(filepath, "neuron_time_F0")
neuron_time_bass = neuron_time_F0;
accuracy_bass = [neuron_time_F0.accuracy]*100;

load(fullfile(base, 'model_comparisons', 'Neuron_Time_F0_Oboe.mat'), ...
	"neuron_time_F0")
accuracy_time = [neuron_time_F0.accuracy]*100;
CFs = [neuron_time_F0.CF];

filename = 'Neuron_Rate_F0_Oboe';
filepath = fullfile(base, 'model_comparisons',filename);
load(filepath, "neuron_rate_F0")
accuracy_rate = [neuron_rate_F0.accuracy]*100;


%% Create figure

figure('Position',[560,563,6.9*ppi,2.6*ppi])
linewidth = 1;
fontsize = 8;
legsize = 7;
labelsize = 12; 

%% A. Confusion matrix of best timing prediction

% Find indices of the best neurons
[temp,originalpos] = sort(accuracy_time, 'descend' );
n = temp(1:3);
best_ind=originalpos(1:3);

% Plot confusion matrix of best neuron
h(1) = subplot(4, 2, 1);
putatives = neuron_time_F0(best_ind(1)).putative;
h(1) = confusionchart(neuron_time_F0(best_ind(1)).C);
%imagesc(neuron_time_F0(best_ind(1)).C);
title(sprintf('Best timing prediction\n%0.02f%%', temp(1)))
set(gca, 'fontsize', fontsize)
%yticks(0:5:10)

%% Second plot

h(2) = subplot(4, 2, 2);
scatter(CFs/1000, accuracy_time, 5, 'filled', 'MarkerEdgeColor','k')
hold on
yline(2.5)
set(gca, 'xscale', 'log')
xlabel('CFs (kHz)')
ylabel('Accuracy (%)')
set(gca, 'fontsize', fontsize)
grid on
xticks([1 2 5 10 20 50 100 200 500])

%% C. Rate vs timing
for itype = 1:2

	h(3+itype) = subplot(4, 2, 3+itype);
	edges = linspace(0, 100, 101);
	chance = 1/length(neuron_rate_F0(1).rate)*100;
	hold on
	if itype == 2
		accuracy = accuracy_time;
		histogram(accuracy,edges, 'Orientation','horizontal')
		yline(chance, 'k', 'LineWidth',linewidth)
		yline(mean(accuracy), 'r', 'LineWidth',linewidth)
		ylim([0 25])
		ylabel('Timing Prediction Accuracy (%)')
	else
		accuracy = accuracy_rate;
		histogram(accuracy,edges)
		xline(chance, 'k', 'LineWidth',linewidth)
		xline(mean(accuracy), 'r', 'LineWidth',linewidth)
		xlim([0 25])
		xlabel('Rate Prediction Accuracy (%)')
	end
	grid on
	set(gca, 'fontsize', fontsize)
end


h(3) = subplot(4, 2, 3);
scatter(accuracy_rate, accuracy_time, 10, 'filled', 'MarkerEdgeColor','k', ...
	'MarkerFaceAlpha',0.5)
hold on
plot([0, 100], [0, 100], 'k')
xlim([0 25])
ylim([0 25])

mdl = fitlm(accuracy_rate, accuracy_time);
x = linspace(0, 100, 20);
y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
plot(x, y, 'r')
yticks(0:5:100)
xticks(0:5:100)
hleg = legend('Neuron', 'Unity', ...
	sprintf('y = %0.2f*x+%0.2f, p=%0.04f', mdl.Coefficients{2, 1}, ...
	mdl.Coefficients{1, 1},mdl.Coefficients{2,4}), 'fontsize', legsize);
hleg.ItemTokenSize = [8, 8];
title('Rate vs Timing Comparison')
set(gca, 'fontsize', fontsize)
xticklabels([])
yticklabels([])
grid on


%% D. Timing: Oboe vs Bassoon

% Find all rows with bassoon and oboe
has_bass = ~cellfun(@isempty, {nat_data.bass_rate});
has_oboe = ~cellfun(@isempty, {nat_data.oboe_rate});
sesh = find(has_bass & has_oboe);
putative = {nat_data(sesh).putative};
nneurons = length(sesh);

for i = 1:nneurons

	putative1 = putative{i};
	ind_bass = find(strcmp(putative1, {neuron_time_bass.putative}));
	ind_oboe = find(strcmp(putative1, {neuron_time_F0.putative}));

	accuracy_oboe(i) = neuron_time_F0(ind_oboe).accuracy*100;
	accuracy_bass1(i) = neuron_time_bass(ind_bass).accuracy*100;
end

for itype = 1:2

	% Plot accuracy of each neuron
	h(6+itype) = subplot(4, 2, 6+itype);
	edges = linspace(0, 100, 31);
	chance = 1/length(neuron_rate_F0(1).rate)*100;
	hold on
	if itype == 2
		accuracy = accuracy_oboe;
		histogram(accuracy,edges, 'Orientation','horizontal')
		yline(chance, 'k', 'LineWidth',linewidth)
		yline(mean(accuracy), 'r', 'LineWidth',linewidth)
		ylim([0 65])
		ylabel('Oboe Prediction Accuracy (%)')
	else
		accuracy = accuracy_bass1;
		histogram(accuracy,edges)
		xline(chance, 'k', 'LineWidth',linewidth)
		xline(mean(accuracy), 'r', 'LineWidth',linewidth)
		xlim([0 65])
		xlabel('Bassoon Prediction Accuracy (%)')
	end
	grid on
	set(gca, 'fontsize', fontsize)
end


h(6) = subplot(4, 2, 6);
scatter(accuracy_bass1, accuracy_oboe, 10, 'filled', 'MarkerEdgeColor','k', ...
	'MarkerFaceAlpha',0.5)
hold on
plot([0, 100], [0, 100], 'k')
xlim([0 65])
ylim([0 65])

mdl = fitlm(accuracy_bass1, accuracy_oboe);
x = linspace(0, 100, 20);
y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
plot(x, y, 'r')
yticks(0:5:100)
xticks(0:5:100)
hleg = legend('Neuron', 'Unity', ...
	sprintf('y = %0.2f*x+%0.2f, p=%0.04f', mdl.Coefficients{2, 1}, ...
	mdl.Coefficients{1, 1},mdl.Coefficients{2,4}), 'fontsize', legsize);
hleg.ItemTokenSize = [8, 8];
title('Oboe vs Bassoon Comparison')
set(gca, 'fontsize', fontsize)
xticklabels([])
yticklabels([])
grid on


%% Position plots
set(h(1), 'Position', [0.07 0.55 0.15 0.32]);
set(h(2), 'Position', [0.07 0.13 0.15 0.28]);


all_fig_positions = ...
	[0.36,0.27,0.24,0.6;...
	0.73,0.27,0.24,0.6]; % left bottom width height

subplot_numbers = [3, 6];
for ipos = 1:2
	fig_position = all_fig_positions(ipos,:);
	nb_position = [fig_position(1),fig_position(2)-0.12,fig_position(3),0.09];
	wb_position = [fig_position(1)-0.04,fig_position(2),0.03,fig_position(4)];
	set(h(subplot_numbers(ipos)), 'Position', fig_position)
	set(h(subplot_numbers(ipos)+1), 'Position', nb_position)
	set(h(subplot_numbers(ipos)+2), 'Position', wb_position)
end

% Create textbox
labels = {'A', 'B', 'C'};
labelleft= [0.01 0.27 0.63];
labelbottom = [repmat(0.96,1, 3) repmat(0.48, 1, 3)];
for ii = 1:3
	annotation('textbox',[labelleft(ii) labelbottom(ii) 0.071 0.058],...
		'String',labels{ii},'FontWeight','bold','FontSize',labelsize,...
		'EdgeColor','none');
end

%% Save figure 

if save_fig == 1
	filename = 'fig5_pitch_single_unit_oboe';
	save_figure(filename)
end
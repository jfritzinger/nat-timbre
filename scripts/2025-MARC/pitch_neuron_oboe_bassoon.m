%% pitch_neuron_oboe_bassoon
clear
save_fig = 1;

%%

[base, ~, ~, ppi] = getPathsNT();
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
accuracy_rate = [neuron_rate_F0.accuracy];

%% Set up figure 

figure('Position',[50 50 4.2*ppi, 5*ppi])
%tiledlayout(3, 4, 'TileIndexing','columnmajor')
linewidth = 2;
fontsize = 18;
titlesize = 20;
legsize = 20;
scattersize = 40;

% D. Timing: Oboe vs Bassoon

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

	accuracy_oboe(i) = neuron_time_F0(ind_oboe).accuracy;
	accuracy_bass1(i) = neuron_time_bass(ind_bass).accuracy;
end

for itype = 1:2

	% Plot accuracy of each neuron
	h(itype+1) = subplot(1, 3, itype+1);
	edges = linspace(0, 1, 31);
	chance = 1/length(neuron_rate_F0(1).rate);
	hold on
	if itype == 2
		accuracy = accuracy_bass1;
		histogram(accuracy,edges, 'Orientation','horizontal')
		yline(chance, 'k', 'LineWidth',linewidth)
		yline(mean(accuracy), 'r', 'LineWidth',linewidth)
		ylim([0 0.65])
		ylabel('Bassoon Timing Accuracy')
	else
		accuracy = accuracy_oboe;
		histogram(accuracy,edges)
		xline(chance, 'k', 'LineWidth',linewidth)
		xline(mean(accuracy), 'r', 'LineWidth',linewidth)
		xlim([0 0.65])
		xlabel('Oboe Timing Accuracy')
		yticks([0 100])
	end
	grid on
	set(gca, 'fontsize', fontsize)
end


h(1) = subplot(1, 3, 1);
scatter(accuracy_oboe, accuracy_bass1, scattersize, 'filled', 'MarkerEdgeColor','k', ...
	'MarkerFaceAlpha',0.5)
hold on
plot([0, 1], [0, 1], 'k', 'LineWidth',2)
xlim([0 0.65])
ylim([0 0.65])

mdl = fitlm(accuracy_oboe,accuracy_bass1);
x = linspace(0, 1, 20);
y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
plot(x, y, '--k', 'LineWidth',2)
yticks(0:0.05:1)
xticks(0:0.05:1)
hleg = legend('Neuron', 'Unity', ...
	sprintf('p=%0.04f', mdl.Coefficients{2,4}), 'fontsize', legsize);
hleg.ItemTokenSize = [20, 20];
title('Timing Comparison')
set(gca, 'fontsize', fontsize)
xticklabels([])
yticklabels([])
grid on


% Position plots
all_fig_positions = [0.32,0.27,0.65,0.65]; % left bottom width height
subplot_numbers = 1;
for ipos = 1
	fig_position = all_fig_positions(ipos,:);
	nb_position = [fig_position(1),fig_position(2)-0.12,fig_position(3),0.09];
	wb_position = [fig_position(1)-0.12,fig_position(2),0.1,fig_position(4)];
	set(h(subplot_numbers(ipos)), 'Position', fig_position)
	set(h(subplot_numbers(ipos)+1), 'Position', nb_position)
	set(h(subplot_numbers(ipos)+2), 'Position', wb_position)
end


% Save figure 
if save_fig == 1
	filename = 'pitch_neuron_oboe_bassoon';
	save_figure_MARC(filename)
end

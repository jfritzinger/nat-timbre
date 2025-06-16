%% pitch_neuron_VS_correlation
clear
save_fig = 1;

%% Load in data 

target = 'Bassoon';
[base, ~, ~, ppi] = getPathsNT();
load(fullfile(base, 'model_comparisons', ['Neuron_Time_F0_' target '.mat']), ...
	"neuron_time_F0")

accuracy_time = [neuron_time_F0.accuracy]*100;
CFs = [neuron_time_F0.CF];
accuracy = [neuron_time_F0.accuracy]*100;

%% Set up figure 

figure('Position',[50,50,4*ppi,4*ppi])
tiledlayout(1, 3)
linewidth = 2;
fontsize = 18;
titlesize = 20;
legsize = 20;
scattersize = 40;
colorsPitch = '#0072B2';

load(fullfile(base,'model_comparisons', 'Data_NT_3.mat'), 'nat_data')

% Get subset of units 
VS_all = zeros(3,1);
for ii = 1:297
	VS_1 = nat_data(ii).bass_VS;
	if ~isempty(VS_1)
		VS_all(ii,:) = mean(VS_1);
	end
end
VS_all(VS_all==0) = [];

scatter(VS_all, accuracy, scattersize, 'filled',...
	'MarkerEdgeColor','k', 'MarkerFaceColor', colorsPitch, 'MarkerFaceAlpha',0.5);
hold on

mdl = fitlm(VS_all, accuracy);
x = linspace(0, 0.6, 20);
y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
plot(x, y, '--k', 'LineWidth',2)
hleg = legend('', ...
	sprintf('p=%0.03f', mdl.Coefficients{2,4}), 'box', 'off', 'Location','northwest');
hleg.ItemTokenSize = [20, 20];
ylabel('Bassoon Timing Accuracy')
xlabel('Vector Strength')
set(gca, 'fontsize', fontsize)
ylim([0 65])
yticks(0:10:60)
grid on

%% Annotate 


%% Save figure 

if save_fig == 1
	filename = 'pitch_neuron_VS_correlation';
	save_figure_MARC(filename)
end

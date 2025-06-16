%% timbre_single_unit_example
clear
save_fig = 0;

%% Load in data

[base, datapath, ~, ppi] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Neuron_Time_Timbre_All.mat'),...
	"neuron_time_timbre")

% Get best unit
accuracy = [neuron_time_timbre.accuracy];
[~,originalpos] = sort(abs(accuracy), 'descend' );
best_ind = originalpos(1:100);
index = best_ind(5);
putative = neuron_time_timbre(index).putative;
%putative = 'R29_TT2_P3_N03';
%index = 149;

% Load in spreadsheet & data
spreadsheet_name = 'Data_Table.xlsx';
sessions = readtable(fullfile(base, spreadsheet_name), ...
	'PreserveVariableNames',true);
load(fullfile(datapath, [putative '.mat']), 'data');

% Find example in spreadsheet
s_ind = strcmp(sessions.Putative_Units, putative);
CF = sessions.CF(s_ind);

%% Set up figure

figure('Position',[50 50 6*ppi, 4.7*ppi])
%tiledlayout(3, 4, 'TileIndexing','columnmajor')
scattersize = 5;
titlesize = 9;
fontsize = 8;
labelsize = 12;
legsize = 7;
linewidth = 1;
spont_color = [0.4 0.4 0.4];
CF_color = [0.7 0.7 0.7];

%% A.

% Plot confusion matrix
h(1) = subplot(3, 4, 1);
C = confusionmat(neuron_time_timbre(index).Instrument, ...
	neuron_time_timbre(index).prediction);
%confusionchart(C)
imagesc(C)
set(gca,'fontsize',fontsize)
xlabel('Prediction Class')
ylabel('True Class')
title('Confusion Matrix', 'fontsize', titlesize)

%% B.
h(2) = subplot(4, 4, 2);

% Load in ERR and mod depth
load(fullfile(base, 'ERR_ModDpth_Bassoon.mat'), "ERR", "modDepth", "pitches")
accurate_predict = diag(C);

% Plot
% scatter(modDepth, accurate_predict, scattersize, 'filled')
% hold on
% mdl = fitlm(modDepth, accurate_predict);
% x = linspace(0, 1, 20);
% y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
% plot(x, y, 'r')
% ylim([0 20])
% hleg = legend('', ...
% 	sprintf('p=%0.04f', mdl.Coefficients{2,4}), 'box', 'off');
% hleg.ItemTokenSize = [8, 8];
xlabel('Modulation Depth')
ylabel('# Correct Classifications')
set(gca,'fontsize',fontsize)
title('Mod Depth', 'fontsize', titlesize)

%% C.

% Load in data
%load(fullfile(base, 'model_comparisons', 'Best_Linear_F0.mat'), "pred_F0", "T")

h(3) = subplot(4, 4, 3);

% r = corrcoef(pred_F0, T.response);
% r2 = r(1, 2)^2;
%
% scatter(T.response, pred_F0, scattersize, 'filled', 'MarkerFaceAlpha',0.5)
% %set(gca, 'xscale', 'log', 'yscale', 'log')
% % ylim([57 588])
% % xlim([57 588])
% hold on
% %title(['R^2 = ' num2str(r2)])
% ylabel('Predicted log10(F0)')
% xlabel('Actual log10(F0)')
% set(gca,'fontsize',fontsize)
% title('Linear Model')


%% D. RM

params_RM = data{2,2}; % Gets binaural response map

% Analysis
data_RM = analyzeRM(params_RM);

% Plot
for ind = 1:3
	h(3+ind) = subplot(4, 4, 3+ind);

	hold on

	if ind == 1
		area(data_RM.freqs,data_RM.rates(:,3),'EdgeColor', '#A49BD0',...
			'FaceColor', '#A49BD0','LineWidth',linewidth) % 33 dB
		spl = '33 dB SPL';
	elseif ind == 2
		area(data_RM.freqs,data_RM.rates(:,4),'EdgeColor', '#5E50A9',...
			'FaceColor', '#5E50A9','LineWidth',linewidth) % 53 dB
		spl = '53 dB SPL';
	elseif ind == 3
		area(data_RM.freqs,data_RM.rates(:,5),'EdgeColor', '#20116B',...
			'FaceColor', '#20116B','LineWidth',linewidth) % 73 dB
		spl = '73 dB SPL';
	end
	plot(data_RM.freqs([1 end]),[1 1]*data_RM.spont,'-','LineWidth',...
		linewidth, 'Color',spont_color)
	xline(CF, '--', 'Color', CF_color,'LineWidth',linewidth);
	text(0.05, 0.90, spl, 'Units', 'normalized', ...
		'VerticalAlignment', 'top', 'FontSize',legsize)
	grid on
	hold off
	ylim([0 max(data_RM.rates, [], 'all')+5])
	yticks([0 40])
	xticks([200 500 1000 2000 5000 10000])
	xticklabels([200 500 1000 2000 5000 10000]./1000)
	set(gca, 'XScale', 'log','fontsize',fontsize);
	if ind == 1
		xlabel('Tone Freq. (kHz)')
		hLegend = legend('', '', 'CF', 'Location','northeast ', ...
			'Box','off');
		hLegend.ItemTokenSize = [12,6];
	elseif ind == 2
		ylabel('Avg. Rate (sp/s)')
		xticklabels([])
	elseif ind == 3
		title('Response Map', 'FontSize', titlesize)
		xticklabels([])
	end
end


%% E. MTF

params_MTF = data{3,2}; % Gets binaural MTFN

% Analysis
data_MTF = analyzeMTF(params_MTF);

% Plot
h(7) = subplot(4, 4, 7);

hold on
line([1 data_MTF.fms(end)], [1 1]*data_MTF.rate(1),'Color',spont_color,...
	'LineWidth',linewidth);
errorbar(data_MTF.fms,data_MTF.rate, ...
	data_MTF.rate_std/sqrt(params_MTF.nrep),'.k', 'LineWidth',...
	linewidth, 'CapSize',4);
line(data_MTF.fms,data_MTF.rate,'Color','k', 'Marker','.', ...
	'MarkerSize',5, 'MarkerFaceColor','w', 'LineWidth', linewidth);
hold off
set(gca, 'XScale', 'log');
xlim([data_MTF.fms(1) data_MTF.fms(end)])
xticks([1 2 5 10 20 50 100 200 500])
xlabel('Modulation Freq. (Hz)')
set(gca,'fontsize',fontsize)
grid on
title('MTF', 'FontSize', titlesize)
ylabel('Avg. Rate (sp/s)')
hLegend = legend('Unmodulated', 'Location','southwest', ...
	'Box','off');
hLegend.ItemTokenSize = [6,6];


%% F. RVF

params_RVF = data{5, 1};

% Plot
h(8) = subplot(4, 4, 8);

if ~isempty(params_RVF)
	data_RVF = analyzeRVF(params_RVF);
	obj = errorbar(data_RVF.velocities,data_RVF.rate,data_RVF.rate_std,'o-', 'LineWidth', 1);
	obj.LineStyle = 'none';
	obj.Marker = 'none';
	obj.Color = 'black';

	nvels = length(data_RVF.velocities);
	max_rvf_rate = max(data_RVF.rate)+data_RVF.rate_std(find(data_RVF.rate==max(data_RVF.rate),1));
	hold on
	plot(data_RVF.velocities(1:nvels/2),data_RVF.rate(1:nvels/2),'k','LineWidth',1)
	plot(data_RVF.velocities(((nvels/2)+1):end),data_RVF.rate(((nvels/2)+1):end),'k','LineWidth',1)
	xline(0,'--', 'LineWidth',2)
	line([1,9],[1 1]*max_rvf_rate*1.04,'Color','blue','LineWidth',1)
	text(4,max_rvf_rate*1.05,'-C','FontSize',9,'Color','blue',...
		'HorizontalAlignment','left','VerticalAlignment','bottom')
	line([-1,-9],[1 1]*max_rvf_rate*1.04,'Color','red','LineWidth',1)
	text(-4,max_rvf_rate*1.05,'+C','FontSize',9,'Color','red',...
		'HorizontalAlignment','right','VerticalAlignment','bottom')
	hold off
	ylabel('Rate (spk/s)')
	xlabel('Velocity (kHz/ms)')
	ylim([0 max_rvf_rate*1.1])
	xlim([min(data_RVF.velocities)*1.05 max(data_RVF.velocities)*1.05])
	xtick = -9:3:9;
	set(gca,'XTick',xtick)
	set(gca,'box','off')
	title('RVF')
	set(gca,'fontsize',fontsize)
end
%% G. PSTH
h(9) = subplot(4, 4, 9);
params_NT = data(7, 2);
data_NT = analyzeNT(params_NT{1});
temporal = analyzeNT_Temporal(data_NT, CF);

% Plot PSTH
max_rate = max(temporal.PSTH, [], 'all');
note_values = round(data_NT.pitch_num);
num_stim = length(note_values);
hold on
for j = 1:num_stim

	% Plot PSTHs
	counts = temporal.PSTH(j,:);
	edges = temporal.t;
	t_bin = edges(1:end-1) + diff(edges)/2; % Bin centers
	x_patch = repelem(edges, 2);
	y_patch = repelem([0; counts(:); 0]', 2);
	y_patch = y_patch(2:end-1); % Creates [0 y1 y1 0 0 y2 y2 0...]
	offset = (j-1)*max_rate; % Adjust offset amount
	patch(x_patch, y_patch + offset, 'b', 'FaceAlpha',0.5, 'EdgeColor','k');

	% Plot smoothed PSTH
	plot(temporal.t(1:end-1), temporal.PSTH_smooth(j,:)+offset,'k', 'LineWidth',1.5);

end
ylim([0 max_rate*num_stim])
xlabel('Time (ms)')
yticks(linspace(max_rate/2, max_rate*num_stim-max_rate/2, num_stim))
yticklabels(note_values)
title('PSTH')
set(gca, 'fontsize', fontsize)


%% H. Period PSTH
h(10) = subplot(4, 4, 10);

% Plot period histogram
max_rate = max(temporal.p_hist, [], 'all');
for j = 1:num_stim

	% Plot PSTHs
	counts = temporal.p_hist(j,:);
	edges = temporal.t_hist(j,:);
	t_bin = edges(1:end-1) + diff(edges)/2; % Bin centers
	x_patch = repelem(edges, 2);
	y_patch = repelem([0; counts(:); 0]', 2);
	y_patch = y_patch(2:end-1); % Creates [0 y1 y1 0 0 y2 y2 0...]
	offset = (j-1)*max_rate; % Adjust offset amount
	patch(x_patch, y_patch + offset, 'b', 'FaceAlpha',0.5, 'EdgeColor','k');
	period = 1/data_NT.pitch_num(j)*1000;
end
ylim([0 max_rate*num_stim])
xlabel('Time (ms)')
yticks(linspace(max_rate/2, max_rate*num_stim-max_rate/2, num_stim))
yticklabels([])
xlim([0 17.5])
%xticks(1:5)
grid on
title('Period Histogram')
set(gca, 'fontsize', fontsize)

%% Arrange plots

left = [0.08 0.35 0.61 0.81];
bottom = linspace(0.08, 0.75, 3);
width = 0.18;
height = 0.19;

set(h(1), 'position', [left(1) bottom(3) width height])
set(h(2), 'position', [left(1) bottom(2) width height])
set(h(3), 'position', [left(1) bottom(1) width height])

set(h(4), 'position', [left(2) bottom(3) width height/3])
set(h(5), 'position', [left(2) bottom(3)+height/3 width height/3])
set(h(6), 'position', [left(2) bottom(3)+height/3*2 width height/3])

set(h(7), 'position', [left(2) bottom(2) width height])
set(h(8), 'position', [left(2) bottom(1) width height])

set(h(9), 'position', [left(3) bottom(1) width 0.87])
set(h(10), 'position', [left(4) bottom(1) width 0.87])

%% Label plots

labelleft= [0 0.3 0.55];
labelbottom = linspace(0.28, 0.94, 3);

annotation('textbox',[labelleft(1) labelbottom(3) 0.071 0.058],...
	'String','A','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(1) labelbottom(2) 0.071 0.058],...
	'String','B','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(1) labelbottom(1) 0.071 0.058],...
	'String','C','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(2) labelbottom(3) 0.071 0.058],...
	'String','D','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(2) labelbottom(2) 0.071 0.058],...
	'String','E','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(2) labelbottom(1) 0.071 0.058],...
	'String','F','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(3) labelbottom(3) 0.071 0.058],...
	'String','G','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');


%% Save figure
%
% if save_fig == 1
% 	filename = 'fig4_2_pitch_single_unit_example';
% 	save_figure(filename)
% end



%% pitch_single_unit_example
clear
save_fig = 0;

%% Load in data

target = 'Bassoon';
[base, datapath, ~, ppi] = getPathsNT();
load(fullfile(base, 'model_comparisons', ['Neuron_Time_F0_' target '.mat']), ...
	"neuron_time_F0")

% Get best unit
accuracy = [neuron_time_F0.accuracy];
[acc_sorted, index] = sort(accuracy, 'descend');

%putative = neuron_time_F0(index(8)).putative;
% putative = 'R29_TT2_P3_N03';
% index = 149;
%putative = 'R29_TT3_P1_N02';
%putative = 'R29_TT3_P2_N11';
%putative = 'R29_TT4_P1_N04';
%putative = 'R29_TT4_P2_N06';
%putative = 'R29_TT3_P3_N02';
%putative = 'R27_TT3_P1_N05';
putative =  'R27_TT2_P11_N03';


% Load in spreadsheet & data
spreadsheet_name = 'Data_Table.xlsx';
sessions = readtable(fullfile(base, spreadsheet_name), ...
	'PreserveVariableNames',true);
load(fullfile(datapath, [putative '.mat']), 'data');

% Find example in spreadsheet
s_ind = strcmp(sessions.Putative_Units, putative);
CF = sessions.CF(s_ind);
F0s = getF0s(target);

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
%rastercolors = {[27,158,119]/256, [217,95,2]/256, [117,112,179]/256};
rastercolors = {[31,88,239]/256, [218,14,114]/256, [255,176,0]/256};

%% D. RM

params_RM = data{2,2}; % Gets binaural response map

% Analysis
data_RM = analyzeRM(params_RM);

% Plot
for ind = 1:3
	h(ind) = subplot(4, 4, ind);

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
h(4) = subplot(4, 4, 4);

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

% params_RVF = data{5, 1};
% data_RVF = analyzeRVF(params_RVF);

% Plot
h(5) = subplot(4, 4, 5);


% obj = errorbar(data_RVF.velocities,data_RVF.rate,data_RVF.rate_std,'o-', 'LineWidth', 1);
% obj.LineStyle = 'none';
% obj.Marker = 'none';
% obj.Color = 'black';
% 
% nvels = length(data_RVF.velocities);
% max_rvf_rate = max(data_RVF.rate)+data_RVF.rate_std(find(data_RVF.rate==max(data_RVF.rate),1));
% hold on
% plot(data_RVF.velocities(1:nvels/2),data_RVF.rate(1:nvels/2),'k','LineWidth',1)
% plot(data_RVF.velocities(((nvels/2)+1):end),data_RVF.rate(((nvels/2)+1):end),'k','LineWidth',1)
% xline(0,'--', 'LineWidth',2)
% line([1,9],[1 1]*max_rvf_rate*1.04,'Color','blue','LineWidth',1)
% text(4,max_rvf_rate*1.05,'-C','FontSize',9,'Color','blue',...
% 	'HorizontalAlignment','left','VerticalAlignment','bottom')
% line([-1,-9],[1 1]*max_rvf_rate*1.04,'Color','red','LineWidth',1)
% text(-4,max_rvf_rate*1.05,'+C','FontSize',9,'Color','red',...
% 	'HorizontalAlignment','right','VerticalAlignment','bottom')
% hold off
% ylabel('Rate (spk/s)')
% xlabel('Velocity (kHz/ms)')
% ylim([0 max_rvf_rate*1.1])
% xlim([min(data_RVF.velocities)*1.05 max(data_RVF.velocities)*1.05])
% xtick = -9:3:9;
% set(gca,'XTick',xtick)
% set(gca,'box','off')
% title('RVF')
% set(gca,'fontsize',fontsize)

%% G. Rasters

h(6) = subplot(4, 4, 6);
hold on
params_NT = data(7, 2);
data_NT = analyzeNT(params_NT{1});
temporal = analyzeNT_Temporal(data_NT, CF);

% Plot dot rasters
num_stim = 40;
for j = 1:num_stim
	offset = (j-1)*20; % Adjust offset amount
	if temporal.VS_p(j)<0.01
		scatter(temporal.x{j}/1000,temporal.y{j}+offset,1, 'filled', 'MarkerFaceColor',rastercolors{mod(j, 3)+1})
	else
		scatter(temporal.x{j}/1000,temporal.y{j}+offset,1, 'filled', 'MarkerFaceColor','k', 'MarkerFaceAlpha',0.5)
	end
	yline(offset, 'k')
end
ylim([0 20*num_stim])
xlabel('Time (ms)')
yticks(linspace(20/2, 20*num_stim-20/2, num_stim))
yticklabels(round(F0s))
title('Rasters')
set(gca,'fontsize',fontsize)
ylabel('F0 (Hz)')
xlim([0 0.15])


%% H. Period PSTH
h(7) = subplot(4, 4, 7);

hold on
max_rate = max([temporal.p_hist{:}]); %max(temporal.p_hist, [], 'all');
for j = 1:num_stim

	% Plot PSTHs
	counts = temporal.p_hist{j,:};
	%counts = smooth_rates(counts,zeros(length(counts), 1),counts+10, 500);
	%counts = zscore(counts);
	edges = temporal.t_hist{j,:};
	t_bin = edges(1:end-1) + diff(edges)/2; % Bin centers
	x_patch = repelem(edges, 2);
	y_patch = repelem([0; counts(:); 0]', 2);
	y_patch = y_patch(2:end-1); % Creates [0 y1 y1 0 0 y2 y2 0...]
	offset = (j-1)*max_rate; % Adjust offset amount
	if temporal.VS_p(j)<0.01
		patch(x_patch, y_patch + offset, rastercolors{mod(j, 3)+1}, 'FaceAlpha',0.8, 'EdgeColor','k');
	else
		patch(x_patch, y_patch + offset, 'k', 'FaceAlpha',0.3, 'EdgeColor','k');
	end
	period = 1/data_NT.F0s_actual(j)*1000;
	yline(offset, 'k')
	plot([period period], [offset j*max_rate], 'k');
end
VS_p = temporal.VS_p;
r_splithalf(:) = temporal.r_splithalf; % Caclulate reliability metric
ylim([0 max_rate*num_stim])
xlabel('Time (ms)')
xticks(0:30)
xticklabels({'0', '', '', '', '', '5', '','','','','10','','','','','15', '',''})
yticks([])
%yticks(linspace(max_rate/2, max_rate*num_stim-max_rate/2, num_stim))
xlim([0 17.5])
yticklabels([])
xtickangle(0)
grid on
title('Period PSTH')
set(gca,'fontsize',fontsize)


%% ISI histograms
h(8) = subplot(4, 4, 8);
[F0s, peak_harm, peak_harm_num] = calcManualF0s(target);

VS_harms2 =flipud(temporal.VS_harms);
peak_harm2 = fliplr(peak_harm_num);
imagesc(1:30, 1:40, VS_harms2)
hold on
for j = 1:num_stim
	rectangle('position', [peak_harm2(j)-0.5 j-0.5, 1, 1], ...
		'EdgeColor','k', 'LineWidth',1)
end
xlim([0.51 12.5])
xlabel('Harmonic Number')
yticklabels([])
c = colorbar;
c.Label.String = 'Vector Strength';
title('        VS to Harmonics')
set(gca,'fontsize',fontsize)

% % Plot ISI histogram
% hold on
% nreps = params_NT{1}.nrep;
% max_rate = max(temporal.ISI_counts_all, [], 'all')-20;
% for j = 1:40
% 	counts = temporal.ISI_counts_all(j,:);
% 	edges = temporal.ISI_edges;
% 
% 	t_bin = edges(1:end-1) + diff(edges)/2; % Bin centers
% 	x_patch = repelem(edges, 2);
% 	y_patch = repelem([0; counts(:); 0]', 2);
% 	y_patch = y_patch(2:end-1); % Creates [0 y1 y1 0 0 y2 y2 0...]
% 	offset = (j-1)*max_rate; % Adjust offset amount
% 	if temporal.VS_p(j)<0.01
% 		patch(x_patch, y_patch + offset,rastercolors{mod(j, 3)+1},  'FaceAlpha',0.8, 'EdgeColor','k');
% 	else
% 		patch(x_patch, y_patch + offset, 'k','FaceAlpha',0.5, 'EdgeColor','k');
% 	end
% 	T = 1/F0s(j)*1000;
% 	plot([T T], [offset (j)*max_rate], 'k');
% end
% ylim([0 max_rate*40])
% xlabel('Time (ms)')
% xticks(0:30)
% xticklabels({'0', '', '', '', '', '5', '','','','','10','','','','','15', '',''})
% xtickangle(0)
% xlim([0 18])
% box on
% yticks([])
% %yticks(linspace(max_rate/2, max_rate*40-max_rate/2, 40))
% yticklabels([])
% grid on
% title('ISI Histogram')
% set(gca,'fontsize',fontsize)
% box off

%% Arrange plots

left = [0.08 0.33 0.53 0.73];
bottom = linspace(0.08, 0.75, 3);
width = 0.16;
height = 0.19;

set(h(1), 'position', [left(1) bottom(3) width height/3])
set(h(2), 'position', [left(1) bottom(3)+height/3 width height/3])
set(h(3), 'position', [left(1) bottom(3)+height/3*2 width height/3])

set(h(4), 'position', [left(1) bottom(2) width height])
set(h(5), 'position', [left(1) bottom(1) width height])

set(h(6), 'position', [left(2) bottom(1) width+0.03 0.87])
set(h(7), 'position', [left(3) bottom(1) width+0.03 0.87])
set(h(8), 'position', [left(4) bottom(1) width+0.02 0.87])

%% Label plots

labelleft= [0 0.27 0.5 0.71];
labelbottom = linspace(0.28, 0.95, 3);

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
annotation('textbox',[labelleft(3) labelbottom(3) 0.071 0.058],...
	'String','E','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(4) labelbottom(3) 0.071 0.058],...
	'String','F','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');


%% Save figure

if save_fig == 1
	filename = 'supp6_pitch_single_unit_example';
	save_figure(filename)
end



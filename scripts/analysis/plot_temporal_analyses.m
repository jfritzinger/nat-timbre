%% plot_temporal_analyses
clear

%% Load in spreadsheet

[base, datapath, ~, ppi] = getPaths();
sheetpath = 'data/2025-manuscript/data-cleaning';
spreadsheet_name = 'PutativeTable2.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);

%% Load in example

% Example session
putative = 'R29_TT3_P5_N02';
filename = sprintf('%s.mat', putative);
load(fullfile(datapath,'neural_data', filename)), 'data';
index = find(cellfun(@(s) strcmp(putative, s), sessions.Putative_Units));
CF = sessions.CF(index);
MTF_shape = sessions.MTF{index};

%% Analysis

% Nat timbre analysis
params_NT = data(14, 2); % 13 is oboe, 14 is bassoon
data_NT = cell(2, 1);
if ~isempty(params_NT{1})
	data_NT = analyzeNT(params_NT{1});
end
temporal = analyzeNT_Temporal(data_NT);

%% Plot PSTH and smoothed PSTH

max_rate = max(temporal.PSTH, [], 'all');
note_values = round(data_NT.pitch_num);

figure('Position',[200,22,403,1233])
tiledlayout(1, 2, 'Padding','compact')

% Plot PSTH full response
nexttile
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
box on
yticks(linspace(max_rate/2, max_rate*num_stim-max_rate/2, num_stim))
yticklabels(note_values)
title('PSTH')
set(gca, 'fontsize', 12)

% Plot period histogram
nexttile
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
box on
yticks(linspace(max_rate/2, max_rate*num_stim-max_rate/2, num_stim))
yticklabels(note_values)
xlim([0 17.5])
%xticks(1:5)
grid on
title('Period Histogram')
set(gca, 'fontsize', 12)

% Plot heatmap
% 
% % Plot as heatmap
% p_hist = temporal.p_hist;
% edges = temporal.t_hist;
% t_bin = edges(1:end-1) + diff(edges)/2; % Bin centers
% t = linspace(0, 5, size(p_hist,2));
% grayMap = [linspace(0, 1, 256)', linspace(0, 1, 256)', linspace(0, 1, 256)'];
% grayMap = flipud(grayMap);
% nexttile
% hh = pcolor(t, data_ST.fpeaks./1000, p_hist);
% set(hh, 'EdgeColor', 'none');
% %colorbar;
% %axis square;
% colormap(grayMap);
% xlim([0 5])
% max_rate = max(p_hist, [], "all");
% clim([0 max_rate])
% ylabel('F0 (Hz)')
% set(gca, 'fontsize', fontsize)

%% Calculate VS sync to 200 Hz

figure
hold on
plot(data_NT.pitch_num, temporal.VS)
%plot(data_NT.pitch_num, smooth_rates(temporal.VS,zeros(num_stim, 1),...
%	ones(num_stim, 1), CF), 'k')
xlabel('Spectral Peak Freq. (Hz)')
ylabel('Vector Strength')
legend('VS')
set(gca, 'fontsize', 12)
title('Example neuron VS, BS')
set(gca, 'xscale', 'log')
xlim([55 600])
xticks([55 110 220 440 600])
ylim([0 1])
grid on
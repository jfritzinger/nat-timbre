%% pitch_single_unit_example 
clear 
save_fig = 1;

%% Load in data

target = 'Bassoon';
[base, datapath, ~, ppi] = getPathsNT();
load(fullfile(base, 'model_comparisons', ['Neuron_Time_F0_' target '.mat']), ...
	"neuron_time_F0")

% Get best unit 
putative = 'R29_TT2_P3_N03';
index = 149;

% Load in spreadsheet & data
spreadsheet_name = 'Data_Table.xlsx';
sessions = readtable(fullfile(base, spreadsheet_name), ...
	'PreserveVariableNames',true);
load(fullfile(datapath, [putative '.mat']), 'data');

% Find example in spreadsheet
s_ind = strcmp(sessions.Putative_Units, putative);
CF = sessions.CF(s_ind);

%% Set up figure 

figure('Position',[50 50 5*ppi, 12*ppi])
%tiledlayout(3, 4, 'TileIndexing','columnmajor')
linewidth = 2;
fontsize = 18;
titlesize = 20;
legsize = 20;
scattersize = 40;
spont_color = [0.4 0.4 0.4];
CF_color = [0.7 0.7 0.7];

%% A. 
% Get bassoon stimulus
tuning = readtable(fullfile(base, 'Tuning.xlsx')); % Load in tuning
target = 'Bassoon';
listing = dir(fullfile(base, 'waveforms', ['*' target '*.wav']));
files = {listing.name};
note_names = extractBetween(files, 'ff.', '.');
[~, index] = ismember(note_names, tuning.Note);
F0s = round(tuning.Frequency(index));
[F0s, ~] = sort(F0s);


%% G. PSTH
h(1) = subplot(1, 2, 1);
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
	patch(x_patch, y_patch + offset, 'b', 'FaceAlpha',1, 'EdgeColor','k');

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
 h(2) = subplot(1, 2, 2);

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
	patch(x_patch, y_patch + offset, 'b', 'FaceAlpha',1, 'EdgeColor','k');
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

left = [0.11 0.57];
bottom = linspace(0.08, 0.75, 3);
width = 0.4;
height = 0.88;

set(h(1), 'position', [left(1) bottom(1) width height])
set(h(2), 'position', [left(2) bottom(1) width height])

%% Annotate 

	annotation('textbox',[0 0.94 0.25 0.058],...
		'String','F0 (Hz)','FontSize',fontsize,...
		'EdgeColor','none');

%% Save figure

if save_fig == 1
	filename = 'pitch_psth_example';
	save_figure_MARC(filename)
end



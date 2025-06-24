%% modelocking_analysis
clear

%% Load in data

target = 'Bassoon';
[base, datapath, ~, ppi] = getPathsNT();
load(fullfile(base, 'model_comparisons', ['Neuron_Time_F0_' target '.mat']), ...
	"neuron_time_F0")




%% Plot
% Get best unit
accuracy = [neuron_time_F0.accuracy];
[ac, ind_high] = sort(accuracy, 'descend');

for ii = 1:10
	putative = neuron_time_F0(ind_high(ii)).putative;
	% putative = 'R29_TT2_P3_N03';
	% index = 149;

	% Load in spreadsheet & data
	spreadsheet_name = 'Data_Table.xlsx';
	sessions = readtable(fullfile(base, spreadsheet_name), ...
		'PreserveVariableNames',true);
	load(fullfile(datapath, [putative '.mat']), 'data');
	% Find example in spreadsheet
	s_ind = strcmp(sessions.Putative_Units, putative);
	CF = sessions.CF(s_ind);
	figure

	% Plot period histogram
	params_NT = data(7, 2);
	data_NT = analyzeNT(params_NT{1});
	temporal = analyzeNT_Temporal(data_NT, CF);
	max_rate = max(temporal.p_hist, [], 'all');
	note_values = round(data_NT.pitch_num);

	num_stim = length(note_values);

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
		hold on

		[pks, locs] = findpeaks(counts/max(temporal.p_hist, [], 'all'), 'MinPeakProminence',0.2);
		scatter(edges(locs), offset+counts(locs), 10, 'filled', 'MarkerFaceColor','k')
		pks_all{j, ii} = pks;
		locs_all{j, ii} = edges(locs);
		num_peaks(j, ii) = length(pks);

	end
	ylim([0 max_rate*num_stim])
	xlabel('Time (ms)')
	yticks(linspace(max_rate/2, max_rate*num_stim-max_rate/2, num_stim))
	yticklabels([])
	xlim([0 17.5])
	%xticks(1:5)
	grid on
	title('Period Histogram')

end

%% 

figure
imagesc(1:40, 1:10, num_peaks')
colorbar

%% Vector strength plots 

load(fullfile(base,'model_comparisons', 'Data_NT_3.mat'), 'nat_data')
[~,best_ind] = sort(abs(accuracy), 'ascend' );

% Get subset of units 
VS_all = zeros(3,40);
for ii = 1:287
	VS_1 = nat_data(best_ind(ii)).bass_VS;
	if ~isempty(VS_1)
		VS_all(ii,:) = VS_1;
	end
end

x = getF0s('Bassoon');
y = 1:287;

pcolor(x, y, VS_all, 'EdgeColor','none', 'EdgeAlpha',0)
set(gca, 'xscale', 'log')
xticks([60 100 200 350 550])
xlabel('F0 (Hz)')
ylabel('All Neurons')
c = colorbar;
title('Vector strength at all neurons',...
	'HorizontalAlignment','center')
c.Label.String = 'Vector strength';
colormap("bone")


%% FUNCTIONS 

function pitch = getF0s(target)

[base, ~, ~, ~] = getPathsNT();
tuning = readtable(fullfile(base, 'Tuning.xlsx')); % Load in tuning
listing = dir(fullfile(base, 'waveforms', '*.wav'));
target_WAV = arrayfun(@(n) contains(listing(n).name, target), 1:numel(listing), 'UniformOutput', false);
wav_nums =  find(cell2mat(target_WAV));
d = dir(fullfile(base,'waveforms', '*.wav'));
all_files = sort({d.name});
nfiles = length(wav_nums);

for i = 1:nfiles
	files{1,i} = all_files{wav_nums(i)};
end

% Sort by frequency of pitch
index = [];
note_names = extractBetween(files, 'ff.','.');
for ii = 1:nfiles % Find index of each note in tuning spreadsheet
	index(ii) = find(strcmp(note_names(ii), tuning.Note));
end
pitch_order = tuning.Frequency(index); % Get freqs of each note
[~, order] = sort(pitch_order); % Sort freqs
pitch = pitch_order(order);

end

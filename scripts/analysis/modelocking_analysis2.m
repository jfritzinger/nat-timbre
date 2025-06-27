%% modelocking_analysis_2

clear

%% Load in data

target = 'Bassoon';
[base, datapath, ~, ppi] = getPathsNT();
load(fullfile(base, 'model_comparisons', ['Neuron_Time_F0_' target '.mat']), ...
	"neuron_time_F0")

% Load in spreadsheet & data
sessions = readtable(fullfile(base, 'Data_Table.xlsx'), ...
	'PreserveVariableNames',true);
F0s = getF0s(target);
accuracy = [neuron_time_F0.accuracy];

%% Plot 10 examples, 5 good and 5 bad

hind = [0, 5];
[acc_sorted, ind_high] = sort(accuracy, 'descend');
r_splithalf = NaN(40, 5);
index = [2, 99, length(ind_high)]; % 18
for ii = 1:3

	figure
	tiledlayout(1, 3, 'Padding','tight', 'TileSpacing','tight')
	putative = neuron_time_F0(ind_high(index(ii))).putative;
	%putative = 'R29_TT2_P3_N03';

	% Load in spreadsheet & data
	load(fullfile(datapath, [putative '.mat']), 'data');
	% Find example in spreadsheet
	s_ind = strcmp(sessions.Putative_Units, putative);
	CF = sessions.CF(s_ind);
	MTF = sessions.MTF(s_ind);
	nexttile
	hold on

	% Analyze data 
	params_NT = data(7, 2);
	data_NT = analyzeNT(params_NT{1});
	temporal = analyzeNT_Temporal(data_NT, CF);

	% Plot dot rasters 
	num_stim = 40;
	for j = 1:num_stim
		offset = (j-1)*20; % Adjust offset amount
		scatter(temporal.x{j}/1000,temporal.y{j}+offset,3, 'filled')
		yline(offset, 'k')
	end
	ylim([0 20*num_stim])
	xlabel('Time (ms)')
	yticks(linspace(20/2, 20*num_stim-20/2, num_stim))
	yticklabels(round(F0s))
	grid on
	title(['Rasters, accuracy = ' num2str(acc_sorted(index(ii)))])


	% Plot period histogram
	%max_rate = max(zscore(temporal.p_hist), [], 'all');
	nexttile
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
			patch(x_patch, y_patch + offset, 'b', 'FaceAlpha',0.5, 'EdgeColor','k');
		else
			patch(x_patch, y_patch + offset, 'k', 'FaceAlpha',0.3, 'EdgeColor','k');
		end
		period = 1/data_NT.pitch_num(j)*1000;
		yline(offset, 'k')
	end
	VS_p = temporal.VS_p;
	r_splithalf(:, ii) = temporal.r_splithalf; % Caclulate reliability metric
	ylim([0 max_rate*num_stim])
	xlabel('Time (ms)')
	xticks(0:30)
	yticks(linspace(max_rate/2, max_rate*num_stim-max_rate/2, num_stim))
	xlim([0 17.5])
	yticklabels(round(VS_p, 3))
	grid on
	title('Period PSTH')

	% Plot ISI histogram 
	nexttile
	hold on
	nreps = params_NT{1}.nrep;
	max_rate = max(temporal.ISI_counts_all, [], 'all');
	for j = 1:40
		counts = temporal.ISI_counts_all(j,:);
		edges = temporal.ISI_edges;

		t_bin = edges(1:end-1) + diff(edges)/2; % Bin centers
		x_patch = repelem(edges, 2);
		y_patch = repelem([0; counts(:); 0]', 2);
		y_patch = y_patch(2:end-1); % Creates [0 y1 y1 0 0 y2 y2 0...]
		offset = (j-1)*max_rate; % Adjust offset amount
		patch(x_patch, y_patch + offset, 'b', 'FaceAlpha',0.5, 'EdgeColor','k');		
		T = 1/F0s(j)*1000;
		plot([T T], [offset (j)*max_rate], 'k');
	end
	ylim([0 max_rate*40])
	xlabel('Time (ms)')
	xticks(0:30)
	xlim([0 18])
	box on
	yticks(linspace(max_rate/2, max_rate*40-max_rate/2, 40))
	yticklabels(round(r_splithalf(:,ii), 2))
	grid on
	title('ISI Histogram')
	set(gca, 'fontsize', 12)

	% Plot spectra with CF marked 

	for j = 1:40
		

	end

end

%% 


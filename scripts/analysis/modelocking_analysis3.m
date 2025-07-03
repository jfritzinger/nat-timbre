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
rastercolors = {[31,88,239]/256, [218,14,114]/256, [255,176,0]/256};

% Get peak harmonic of each stimulus 
[F0s, peak_harm] = calcManualF0s(target);

%% Plot 10 examples, 5 good and 5 bad

hind = [0, 5];
[acc_sorted, ind_high] = sort(accuracy, 'descend');
r_splithalf = NaN(40, 5);
index = 18; %[2, 99, length(ind_high)]; % 18
for ii = 1 %:3

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
		if temporal.VS_p(j)<0.01
			scatter(temporal.x{j}/1000,temporal.y{j}+offset,3, 'filled', 'MarkerFaceColor',rastercolors{mod(j, 3)+1})
		else
			scatter(temporal.x{j}/1000,temporal.y{j}+offset,3, 'filled', 'MarkerFaceColor','k', 'MarkerFaceAlpha',0.9)
		end
		yline(offset, 'k')
	end
	ylim([0 20*num_stim])
	xlabel('Time (ms)')
	yticks(linspace(20/2, 20*num_stim-20/2, num_stim))
	yticklabels(round(F0s))
	grid on
	xlim([0 0.15])
	title(['Rasters, accuracy = ' num2str(acc_sorted(index(ii)))])


	% Plot period histogram
	%max_rate = max(zscore(temporal.p_hist), [], 'all');
	nexttile
	hold on
	max_rate = max([temporal.p_hist{:}])-max([temporal.p_hist{:}])/3; %-150; %max(temporal.p_hist, [], 'all');
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
		%if temporal.VS_p_harmns(j, 3)<0.01
			patch(x_patch, y_patch + offset, rastercolors{mod(j, 3)+1}, 'FaceAlpha',0.9, 'EdgeColor','k');
		else
			patch(x_patch, y_patch + offset, 'k', 'FaceAlpha',0.9, 'EdgeColor','k');
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
	%xlim([0 17.5])
	yticklabels(round(VS_p, 3))
	%yticklabels(round(temporal.VS_p_harmns(:, 2), 3))
	grid on
	title('Period PSTH')

	% figure
	% hold on
	% for j = 1:num_stim
	% 	VS_harms = temporal.VS_harms(j, :);
	% 	p_VS_harms = temporal.VS_p_harmns(j,:);
	% 	sig = p_VS_harms<0.01;
	% 	yline(j, 'k')
	% 	%plot(VS_harms+j)
	% 	for jj = 1:10
	% 		if sig(jj)==1
	% 			scatter(jj, VS_harms(jj)+j, 8, 'k', 'filled')
	% 		end
	% 	end
	% 
	% 	x_patch = repelem(0:11, 2);
	% 	y_patch = repelem([0 VS_harms 0], 2);
	% 	patch(x_patch, y_patch + j, rastercolors{mod(j, 3)+1}, 'FaceAlpha',0.9, 'EdgeColor','k');
	% end
	% xlim([1 10])
	% ylim([1 41])
	% yticklabels([])
	% xlabel('Harmonic #')
	% title('VS to each harmonic')

	% Get vector strength for each harmonic closest to CF 
	figure
	tiledlayout(5, 8)
	for j = 1:num_stim
		nexttile
		harms = temporal.harms(j,:);
		VS_harms = temporal.VS_harms(j,:);
		[~, near_ind] = min(abs(CF-harms));
		nearest_harm = harms(near_ind);
		stem(harms, VS_harms);
		hold on
		xline(CF)
		xline(peak_harm(j), 'r')
		title(['F0=' num2str(round(harms(1)))])
		ylim([0 0.8])		
	end

	% Create image plot of phase locking at different harmonics 
	figure
	tiledlayout(40, 1, "TileSpacing","none")
	harms2 = flipud(temporal.harms);
	VS_harms2 =flipud(temporal.VS_harms);
	for j = 1:num_stim
		nexttile
		imagesc(harms2(j,:), 1, VS_harms2(j,:))
		xline(peak_harm(j), 'r', 'LineWidth',2)
		xlim([0 2500])
		clim([0 0.8])
	end

	% Plot ISI histogram 
	nexttile
	hold on
	nreps = params_NT{1}.nrep;
	max_rate = max(temporal.ISI_counts_all, [], 'all'); %-50;
	for j = 1:40
		counts = temporal.ISI_counts_all(j,:);
		edges = temporal.ISI_edges;

		t_bin = edges(1:end-1) + diff(edges)/2; % Bin centers
		x_patch = repelem(edges, 2);
		y_patch = repelem([0; counts(:); 0]', 2);
		y_patch = y_patch(2:end-1); % Creates [0 y1 y1 0 0 y2 y2 0...]
		offset = (j-1)*max_rate; % Adjust offset amount
		if temporal.VS_p(j)<0.01
			patch(x_patch, y_patch + offset, rastercolors{mod(j, 3)+1}, 'FaceAlpha',0.9, 'EdgeColor','k');
		else
			patch(x_patch, y_patch + offset, 'k', 'FaceAlpha',0.9, 'EdgeColor','k');
		end
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


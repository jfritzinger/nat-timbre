%% modelocking_analysis
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

figure
tiledlayout(1, 10, 'Padding','tight', 'TileSpacing','none')
hind = [0, 5];
for igood = 1
	if igood == 1
		[~, ind_high] = sort(accuracy, 'descend');
	else
		[~, ind_high] = sort(accuracy, 'ascend');
	end
	r_splithalf = NaN(40, 5);
	for ii = 1:10
		putative = neuron_time_F0(ind_high(87+ii)).putative;

		% Load in spreadsheet & data
		load(fullfile(datapath, [putative '.mat']), 'data');
		% Find example in spreadsheet
		s_ind = strcmp(sessions.Putative_Units, putative);
		CF = sessions.CF(s_ind);
		nexttile

		% Plot period histogram
		params_NT = data(7, 2);
		data_NT = analyzeNT(params_NT{1});
		temporal = analyzeNT_Temporal(data_NT, CF);
		%max_rate = max(zscore(temporal.p_hist), [], 'all');
		max_rate = max(temporal.p_hist, [], 'all');
		note_values = round(data_NT.pitch_num);
		num_stim = length(note_values);
		for j = 1:num_stim

			% Plot PSTHs
			counts = temporal.p_hist(j,:);
			%counts = smooth_rates(counts,zeros(length(counts), 1),counts+10, 500);
			%counts = zscore(counts);

			edges = temporal.t_hist(j,:);
			t_bin = edges(1:end-1) + diff(edges)/2; % Bin centers
			x_patch = repelem(edges, 2);
			y_patch = repelem([0; counts(:); 0]', 2);
			y_patch = y_patch(2:end-1); % Creates [0 y1 y1 0 0 y2 y2 0...]
			offset = (j-1)*max_rate; % Adjust offset amount
			patch(x_patch, y_patch + offset, 'b', 'FaceAlpha',0.5, 'EdgeColor','k');
			period = 1/data_NT.pitch_num(j)*1000;
			hold on

			% Caclulate reliability metric 
			r_splithalf(j, ii) = temporal.r_splithalf(j);

			if r_splithalf(j, ii) >= 0.2
				[pks, locs] = findpeaks(counts/max(temporal.p_hist, [], 'all'), 'MinPeakProminence',0.1);
				scatter(edges(locs), offset+counts(locs), 5, 'filled', 'MarkerFaceColor','r')
				pks_all{j, ii} = pks;
				locs_all{j, ii} = edges(locs);
				num_peaks(j, ii) = length(pks);
			else
				pks_all{j, ii} = [];
				locs_all{j, ii} = [];
				num_peaks(j, ii) = 0; 
			end
		end
		ylim([0 max_rate*num_stim])
		xlabel('Time (ms)')
		xticks(0:30)
		yticks(linspace(max_rate/2, max_rate*num_stim-max_rate/2, num_stim))
		yticklabels(round(r_splithalf(:,ii), 2))
		%xlim([0 17.5])
		%xticks(1:5)
		grid on

	end
end

%% Plot number of peaks from peak picking

% Get best unit
accuracy = [neuron_time_F0.accuracy];
[ac, ind_high] = sort(accuracy, 'ascend');

for ii = 1:length(accuracy)

	putative = neuron_time_F0(ind_high(ii)).putative;
	load(fullfile(datapath, [putative '.mat']), 'data');
	s_ind = strcmp(sessions.Putative_Units, putative);
	CF = sessions.CF(s_ind);

	% Plot period histogram
	params_NT = data(7, 2);
	data_NT = analyzeNT(params_NT{1});
	temporal = analyzeNT_Temporal(data_NT, CF);
	max_rate = max(temporal.p_hist, [], 'all');
	note_values = round(data_NT.pitch_num);
	num_stim = length(note_values);
	for j = 1:num_stim
		counts = temporal.p_hist(j,:);
		edges = temporal.t_hist(j,:);
		period = 1/data_NT.pitch_num(j)*1000;

		% Caclulate reliability metric
		r_splithalf(j, ii) = temporal.r_splithalf(j);
		if r_splithalf(j, ii) >= 0.2
			[pks, locs] = findpeaks(counts/max(temporal.p_hist, [], 'all'), 'MinPeakProminence',0.1);
			pks_all{j, ii} = pks;
			locs_all{j, ii} = edges(locs);
			num_peaks(j, ii) = length(pks);
		else
			pks_all{j, ii} = [];
			locs_all{j, ii} = [];
			num_peaks(j, ii) = 0;
		end

		% VS
		VS_all2(j, ii) = temporal.VS(j);
	end
end

%% Mode locking plots 

figure
pcolor(F0s, 1:length(accuracy), num_peaks', 'EdgeColor','none', 'EdgeAlpha',0)
colorbar
set(gca, 'xscale', 'log')

figure
mean_peaks = mean(num_peaks, 1);
scatter(1:length(accuracy), mean_peaks, 'filled')
ylabel('Mean # Peaks')
xlabel('Neurons, low -> high accuracy')

figure
mean_peaks = mean(num_peaks, 1);
scatter(ac, mean_peaks, 'filled')
ylabel('Mean # Peaks')
xlabel('Accuracy')

%% Reliability plots 

figure
nexttile
pcolor(F0s, 1:length(accuracy), r_splithalf', 'EdgeColor','none', 'EdgeAlpha',0)
c = colorbar;
set(gca, 'xscale', 'log')
xlabel('F0s (Hz)')
ylabel('Neurons, low -> high accuracy')
title('Reliability')
c.Label.String = 'Reliability';

nexttile
mean_splithalf = mean(r_splithalf, 1);
scatter(1:length(accuracy), mean_splithalf, 'filled')
ylabel('Reliability')
xlabel('Neurons, low -> high accuracy')

nexttile
mean_splithalf = mean(r_splithalf, 1);
scatter(ac, mean_splithalf, 15, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.5)
ylabel('Average Reliability')
xlabel('Accuracy')
title('Avg Reliability vs F0 accuracy (bassoon)')

%% Vector strength plots 

figure
nexttile
pcolor(F0s, 1:length(accuracy), VS_all2', 'EdgeColor','none', 'EdgeAlpha',0)
colorbar
set(gca, 'xscale', 'log')
xlabel('F0s (Hz)')
ylabel('Neurons, low -> high accuracy')
title('Vector Strength')

nexttile
mean_VS = mean(VS_all2, 1);
scatter(1:length(accuracy), mean_VS, 'filled')
ylabel('Vector Strength')
xlabel('Neurons, low -> high accuracy')

nexttile
mean_VS = mean(VS_all2, 1);
scatter(ac, mean_VS, 15, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.5)
ylabel('Average Vector Strength')
xlabel('Accuracy')
title('Avg VS vs F0 accuracy (bassoon)')

%% Plot correlations as histograms 

[~,best_ind] = sort(abs(accuracy), 'ascend' );
for ii = 1:287
	C_acc(ii,:) = diag(neuron_time_F0(best_ind(ii)).C);
end

for ii = 1:287

	VS_C = corrcoef(C_acc(ii,:), VS_all2(:,ii));
	r_VS_C(ii) = VS_C(1, 2);
	rel_C = corrcoef(C_acc(ii,:), r_splithalf(:,ii));
	r_rel_C(ii) = rel_C(1, 2);

end

figure
boxplot([r_VS_C; r_rel_C]')
hold on
swarmchart(ones(287, 1), r_VS_C)
swarmchart(ones(287, 1)*2, r_rel_C)

r_c = [r_VS_C; r_rel_C]';
[h,p,ci,stats] = ttest(r_c);

ylabel('Correlation between __ and accuracy')
xlabel('Groups')
xticklabels({'VS', 'Reliability'})

%% Plot PSTH for good example

figure()
tiledlayout(1, 2, 'TileSpacing','none')
h(1) = nexttile;
for igood = 1
	if igood == 1
		[~, ind_high] = sort(accuracy, 'descend');
	else
		[~, ind_high] = sort(accuracy, 'ascend');
	end
	r_splithalf = NaN(40, 5);
	for ii = 2
		putative = neuron_time_F0(ind_high(ii)).putative;

		% Load in spreadsheet & data
		load(fullfile(datapath, [putative '.mat']), 'data');
		% Find example in spreadsheet
		s_ind = strcmp(sessions.Putative_Units, putative);
		CF = sessions.CF(s_ind);

		% Plot period histogram
		params_NT = data(7, 2);
		data_NT = analyzeNT(params_NT{1});
		temporal = analyzeNT_Temporal(data_NT, CF);
		%max_rate = max(zscore(temporal.p_hist), [], 'all');
		max_rate = max(temporal.p_hist, [], 'all');
		note_values = round(data_NT.pitch_num);
		num_stim = length(note_values);
		for j = 1:num_stim

			% Plot PSTHs
			counts = temporal.p_hist(j,:);
			%counts = smooth_rates(counts,zeros(length(counts), 1),counts+10, 500);
			%counts = zscore(counts);

			edges = temporal.t_hist(j,:);
			t_bin = edges(1:end-1) + diff(edges)/2; % Bin centers
			x_patch = repelem(edges, 2);
			y_patch = repelem([0; counts(:); 0]', 2);
			y_patch = y_patch(2:end-1); % Creates [0 y1 y1 0 0 y2 y2 0...]
			offset = (j-1)*max_rate; % Adjust offset amount
			patch(x_patch, y_patch + offset, 'b', 'FaceAlpha',0.5, 'EdgeColor','k');
			period = 1/data_NT.pitch_num(j)*1000;
			hold on

			% Caclulate reliability metric 
			r_splithalf(j, ii) = temporal.r_splithalf(j);

			if r_splithalf(j, ii) >= 0.2
				[pks, locs] = findpeaks(counts/max(temporal.p_hist, [], 'all'), 'MinPeakProminence',0.1);
				scatter(edges(locs), offset+counts(locs), 15, 'filled', 'MarkerFaceColor','r')
				pks_all{j, ii} = pks;
				locs_all{j, ii} = edges(locs);
				num_peaks(j, ii) = length(pks);
			else
				pks_all{j, ii} = [];
				locs_all{j, ii} = [];
				num_peaks(j, ii) = 0; 
			end
		end
		ylim([0 max_rate*num_stim])
		xticks(0:30)
		xlabel('Time (ms)')
		yticks(linspace(max_rate/2, max_rate*num_stim-max_rate/2, num_stim))
		yticklabels(F0s)
		%xlim([0 17.5])
		%xticks(1:5)
		grid on
	end
end
title('Period PSTH')

%% Calculate ISI for each  for a good example 
figure

nreps = params_NT{1}.nrep;
ISI_all = cell(1, 40);
nbins = 80;
edges = linspace(0, 20, 81);
counts_all = zeros(40, nbins); % Pre-allocate counts_all
for jj = 1:40
	x = temporal.x{jj} / 1000; % ms % spike times 
	y = temporal.y{jj}; % Reps 

	ISI = arrayfun(@(ii) diff(x(y == ii)), 1:nreps, 'UniformOutput', false);
	ISI_all{jj} = vertcat(ISI{:});
	counts_all(jj, :) = histcounts(ISI_all{jj}, edges);
end

% Plot ISI histograms
h(2) = nexttile;
max_rate = max(counts_all, [], 'all')-50;
for j = 1:40
	counts = counts_all(j,:);

	t_bin = edges(1:end-1) + diff(edges)/2; % Bin centers
	x_patch = repelem(edges, 2);
	y_patch = repelem([0; counts(:); 0]', 2);
	y_patch = y_patch(2:end-1); % Creates [0 y1 y1 0 0 y2 y2 0...]
	offset = (j-1)*max_rate; % Adjust offset amount
	patch(x_patch, y_patch + offset, 'b', 'FaceAlpha',0.5, 'EdgeColor','k');
	hold on
	T = 1/F0s(j)*1000;
	plot([T T], [offset (j)*max_rate], 'k');
end
ylim([0 max_rate*40])
xlabel('Time (ms)')
xticks(0:30)
box on
yticks(linspace(max_rate/2, max_rate*40-max_rate/2, 40))
yticklabels(round(r_splithalf(:,ii), 2))
grid on
title('ISI Histogram')
set(gca, 'fontsize', 12)

%% Calculate ISI scatter plot

% Calculate ISI for each rep
ISI_all = cell(1, 40);
counts_all = zeros(40, 50); % Pre-allocate counts_all
for jj = 1:40
	ISI = [];
	x_isi = [];
	y_isi = [];
	x = temporal.x{jj}/1000; % ms
	y = temporal.y{jj};
	for ii = 1:20 % for each repetition
		valid = y==ii;
		isi = [0; diff(x(valid))];
		ISI = [ISI isi'];
		x_isi = [x_isi isi(1:end-1)'];
		y_isi = [y_isi isi(2:end)'];
		index = find(valid);
		rep_ind(ii) = index(1);
	end
	rep_ind_all{jj} = rep_ind;
	ISI_all{jj} = ISI;
	x_isi_all{jj} = x_isi;
	y_isi_all{jj} = y_isi;
end


figure('Position',[560,42,1050,806])
tiledlayout(6, 7, 'Padding','compact')
for jj = 1:40
	T = 1/F0s(jj)*1000;
	nexttile
	hold on
	line(repmat([0;10*T],1,10),[1;1]*(1:10)*T,'Color','r', 'linewidth', 0.3)
	line([1;1]*(1:10)*T,repmat([0;10*T],1,10),'Color','r', 'linewidth', 0.3)
	axis square
	scatter(x_isi_all{jj}, y_isi_all{jj},5, 'filled', 'k', 'markerfacealpha', 0.5)

	if ismember(jj, 37:41)
		xlabel('ISI n (ms)')
	end

	if ismember(jj, 1:7:41)
		ylabel('ISI n+1 (ms)')
	end
	title(sprintf('%0.0f Hz', F0s(jj)))
	xlim([0 5*T])
	ylim([0 5*T])
end

% %% A. Significance Test 1
% % To assess whether ISI dependencies are enough to produce phase-locking, 
% % we tested the surrogate spike train using a standard Rayleigh test. 
% % We denote by Ris the value of this test in the rest of this article.
% 
% for ii = 1:40
% 
% 	% An Ris score of 13.8 (Rayleigh criterion) and a z score of 2 is
% 	% required to meet our criterion for "mode-locking."
% 	original_isi = ISI_all{ii};
% 	period = 1/F0s(ii)*1000;
% 
% 	spike_times = temporal.x{ii}/1000;
% 	rep_ind = [rep_ind_all{ii} length(original_isi)];
% 	spike_times3 = [];
% 	for irep = 1:20
% 		spike_times2 = cumsum([spike_times(rep_ind(irep)), ...
% 			original_isi(rep_ind(irep)+1:rep_ind(irep+1)-1)])';
% 		spike_times3 = [spike_times3; spike_times2];
% 	end
% 
% 	phases = 2 * pi * mod(spike_times, period) / period;
% 	phases2 = 2 * pi * mod(spike_times3, period) / period;
% 
% 	% figure
% 	% nexttile
% 	% plot(spike_times)
% 	% hold on
% 	% plot(spike_times3)
% 	% 
% 	% nexttile
% 	% histogram(phases, 31)
% 	% hold on
% 	% histogram(phases2', 31)
% 
% 	% Create surrogate ISI sequence
% 	surrogate_isi = create_surrogate_isi(original_isi, rep_ind);
% 
% 	% Plot ISI histogram and shuffled histogram
% 	% figure
% 	% edges = linspace(0, 20, 81);
% 	% histogram(original_isi, edges, 'FaceAlpha',0.5);
% 	% hold on
% 	% histogram(surrogate_isi, edges, 'FaceAlpha',0.5);
% 
% 	% Perform Rayleigh test
% 	surrogate_spike_times = cumsum([0, surrogate_isi]);
% 	surrogate_phases = 2 * pi * mod(surrogate_spike_times, period) / period;
% 	[p_value(ii), Ris(ii)] = circ_rtest(phases2);
% end
% passed = sum(p_value<0.05);
% 
% % But none have Ris > 13.8
% 
% %% B. Significance Test 2
% 
% % To test condition B, a surrogate spike train is created by randomly 
% % shuffling the phases of the genuine spike train (see Fig. 3B). 
% 
% % Shuffle the phases of the spike train
% 
% % Plot the period histogram and shuffled period histogram
% 
% 
% % The dissimilarity of ISI sequences for surrogate versus genuine spike 
% % trains was assessed by computing the root mean squared error (RMSe)
% % between their ISI distributions. We compared this RMSe from the 
% % phase shuffling to the RMSe within data trials. 
% % RMSEps = rmse();
% % RMSEtrials = rmse();
% 
% % The resulting statistic is expressed in the results as a z score (Zps)
% % Zps = (RMSEps - RMSEtrials) / std(trials);
% 
% 
% % Check if z-score is greater than 2
% % if z_score > 2
% % 	disp('Z-score passed');
% % else
% % 	disp('Z-score did not pass');
% % end
% 
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function surrogate_isi = create_surrogate_isi(original_isi, rep_ind)
% % 1) randomly pick from the data one ISI, say ISIn, and assign ISIn to the 
% % surrogate sequence. 
% % 2) Find {kISIk  ISIn}, i.e., all the ISIs equal to 
% % ISIn, and build the list of ISI pairs: L1  {(ISIk, ISIk 1)}. And 
% % 3) randomly pick one pair from L1 , say (ISIm , ISIm  1 ) and 
% % assign ISIm 1 to the surrogate sequence.
% % Then we repeat steps 1–3 until no more pairs remain. So we find in the 
% % data all ISIs that are equal to ISIm 1 build L2  {(ISIi, ISIi 1)}, 
% % choose a pair in L2 and assign ISIi 1 to surrogate sequence. 
% 
% 
% % surrogate_isi = [];
% % remaining_isi = original_isi;
% % 
% % while ~isempty(remaining_isi)
% % 
% % 	% Step 1: Randomly pick an ISI
% % 	idx = randi(length(remaining_isi));
% % 	isi_n = remaining_isi(idx);
% % 	surrogate_isi = [surrogate_isi, isi_n];
% % 	remaining_isi(idx) = [];
% % 
% % 	% Step 2: Find all ISIs equal to ISI_n and build pairs
% % 	equal_indices = find(remaining_isi == isi_n);
% % 	if ~isempty(equal_indices)
% % 		pairs = [remaining_isi(equal_indices), remaining_isi(equal_indices + 1)];
% % 
% % 		% Step 3: Randomly pick a pair and assign ISI_m+1
% % 		pair_idx = randi(size(pairs, 1));
% % 		isi_m_plus_1 = pairs(pair_idx, 2);
% % 		surrogate_isi = [surrogate_isi, isi_m_plus_1];
% % 
% % 		% Remove used ISIs from remaining_isi
% % 		remaining_isi(equal_indices(pair_idx)) = [];
% % 		remaining_isi(equal_indices(pair_idx)) = [];
% % 	end
% % end
% 
% % Step 1: Initialize surrogate sequence and select random starting ISI
% surrogate_isi = zeros(size(original_isi));
% start_idx = randi(length(original_isi));
% surrogate_isi(1) = original_isi(start_idx);
% current_isi = surrogate_isi(1);
% 
% % Precompute all consecutive ISI pairs from original data
% all_pairs = [original_isi(1:end-1)', original_isi(2:end)'];
% 
% for i = 2:length(original_isi)
% 	% Step 2: Find all pairs where first element equals current ISI
% 	valid_pairs = all_pairs(all_pairs(:,1) == current_isi, :);
% 
% 	if isempty(valid_pairs)
% 		% If no pairs found, select random ISI from original data
% 		current_isi = original_isi(randi(length(original_isi)));
% 	else
% 		% Step 3: Randomly select a pair and take its next ISI
% 		selected_pair = valid_pairs(randi(size(valid_pairs,1)), :);
% 		current_isi = selected_pair(2);
% 	end
% 
% 	% Assign to surrogate sequence
% 	surrogate_isi(i) = current_isi;
% end
% end
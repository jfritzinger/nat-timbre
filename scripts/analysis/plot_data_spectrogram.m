% plot_rate_spectrogram
clear

%% Load in spreadsheet
[base, datapath, savepath, ppi] = getPaths();
modelpath = '/Volumes/Nat-Timbre/data/manuscript';
sheetpath = 'scripts/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);

%% Get a list of sessions

nat(:,1) = cellfun(@(s) contains(s, 'R'), sessions.Bassoon);
nat(:,2) = cellfun(@(s) contains(s, 'R'), sessions.Oboe);
all_sets = nat(:,1) & nat(:,2);
table = sessions(all_sets, :);

%% Load in each dataset and plot rate response

colors = {"#0072BD", "#D95319"};
j = 4;
for j = 1

	% Load in data
	%putative = table.Putative_Units{j};
	%CF = table.CF(j);
	MTF_shape = table.MTF{j};

	putative = 'R25_TT3_P8_N08'; % or 08
	CF = 6063;
	%putative = 'R25_TT4_P8_N15';
	%CF = 9190;
	load(fullfile(datapath, 'neural_data', [putative '.mat']))

	% Analyze data
	param = data(13:14, 2);
	data_NT = cell(2, 1);
	for ii = 1
		data_NT{ii} = analyzeNT(param{ii});
	end
	overlap = data_NT{2}.pitch_num(25:end);

	% Plot
	% figure('Position',[130,520,1121,273])
	% hold on
	% for ii = 1:16
	% 	ind1 = data_NT{1}.pitch_num == overlap(ii);
	% 	ind2 = data_NT{2}.pitch_num == overlap(ii);
	% 	plot([overlap(ii) overlap(ii)], [data_NT{1}.rate(ind1) data_NT{2}.rate(ind2)], 'k')
	% end
	% scatter(data_NT{1}.pitch_num, data_NT{1}.rate, 'filled', 'MarkerFaceColor',colors{1}, 'MarkerEdgeColor','k');
	% scatter(data_NT{2}.pitch_num, data_NT{2}.rate, 'filled', 'MarkerFaceColor',colors{2}, 'MarkerEdgeColor','k');
	% plot(data_NT{1}.pitch_num(1:17), data_NT{1}.rate(1:17), 'Color',colors{1})
	% plot(data_NT{2}.pitch_num(25:end), data_NT{2}.rate(25:end), 'Color',colors{2 })
	% 
	% set(gca, 'xscale', 'log')
	% xticks([58 116 233 466 932])
	% grid on
	% title(sprintf('%s MTF, CF = %0.0f Hz', MTF_shape, CF))
	% xlabel('F0 (Hz)')
	% ylabel('Spike Rate (sp/s')
	% set(gca, 'FontSize', 14)

	%%

	jj = 1;
	[Timbres,~,Timbrei] = unique([param{jj}.list.iwav_file].');
	num_Timbres = length(Timbres);

	t_spike_rel = data_NT{jj}.t_spike_rel;
	rel_id = data_NT{jj}.rel_id;
	rep = data_NT{jj}.rep;

	figure;
	tiledlayout(num_Timbres, 1, 'TileSpacing','none')
	for ii = 1:num_Timbres

		j = data_NT{jj}.j(ii,:);
		x = t_spike_rel(ismember(rel_id,j));
		x = x/1000000;
		ind2 = x>0.3;
		x(ind2) = [];
		y = rep(rel_id(ismember(rel_id,j)));
		y(ind2) = [];

		for i = 1:20
			ind = y == i;
			spikeTimes{i} = x(ind);
		end

		% Define parameters
		F0 = data_NT{jj}.pitch_num(ii);
		periodLength = 1/F0; % Length of one period in seconds
		numPeriods = 5; % Number of periods to show
		binSize = (1/F0)/25; % Bin size in seconds
		numRepetitions = 20;

		% Create time bins
		timeBins = 0:binSize:(periodLength*numPeriods);

		% Count spikes in each bin across all repetitions
		spikeCount = zeros(size(timeBins));
		for i = 1:numRepetitions
			trialSpikes = spikeTimes{i};
			for j = 1:length(trialSpikes)

				% Wrap spike times to fit within the specified number of periods
				wrappedTime = mod(trialSpikes(j), periodLength*numPeriods);
				binIndex = find(timeBins <= wrappedTime, 1, 'last');
				spikeCount(binIndex) = spikeCount(binIndex) + 1;
			end
		end

		% Calculate firing rate (spikes per second)
		firingRate = spikeCount / (binSize * numRepetitions);

		% Plot the PSTH
		nexttile
		scatter(x, y, '.k')
		%bar(timeBins, spikeCount, 'histc');
		%xlim([0 periodLength*numPeriods]);
		yticks([])

		% Add vertical lines to separate periods
		for i = 1:numPeriods-1
			xline(i*periodLength, '--r', 'LineWidth',1.5);
		end

		% Create matrix
		t = timeBins;
		spec(ii,:) = firingRate;
	end

	harm_num = CF./data_NT{jj}.pitch_num;
	figure
	pcolor(t, harm_num, spec);
	shading interp
	colorbar

end
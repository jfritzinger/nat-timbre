%% plot_overlapping_results.m
%
% Script that plots rate responses to oboe and bassoon from the same neuron
% and overlaps the responses that have the same F0s
%
% Author: J. Fritzinger
% Created: 2022-10-14; Last revision: 2024-10-14
%
% -------------------------------------------------------------------------
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
for j = 1:20

	% Load in data
	putative = table.Putative_Units{j};
	CF = table.CF(j);
	MTF_shape = table.MTF{j};
	load(fullfile(datapath, 'neural_data', [putative '.mat']))

	% Analyze data
	param = data(13:14, 2);
	data_NT = cell(2, 1);
	for ii = 1:2
		data_NT{ii} = analyzeNT(param{ii});
	end
	overlap = data_NT{2}.pitch_num(25:end);


	% Plot
	figure('Position',[130,520,1121,273])
	hold on
	for ii = 1:16
		ind1 = data_NT{1}.pitch_num == overlap(ii);
		ind2 = data_NT{2}.pitch_num == overlap(ii);
		plot([overlap(ii) overlap(ii)], [data_NT{1}.rate(ind1) data_NT{2}.rate(ind2)], 'k')
	end
	scatter(data_NT{1}.pitch_num, data_NT{1}.rate, 'filled', 'MarkerFaceColor',colors{1}, 'MarkerEdgeColor','k');
	scatter(data_NT{2}.pitch_num, data_NT{2}.rate, 'filled', 'MarkerFaceColor',colors{2}, 'MarkerEdgeColor','k');
	plot(data_NT{1}.pitch_num(1:17), data_NT{1}.rate(1:17), 'Color',colors{1})
	plot(data_NT{2}.pitch_num(25:end), data_NT{2}.rate(25:end), 'Color',colors{2 })

	set(gca, 'xscale', 'log')
	xticks([58 116 233 466 932])
	grid on
	title(sprintf('%s MTF, CF = %0.0f Hz', MTF_shape, CF))
	xlabel('F0 (Hz)')
	ylabel('Spike Rate (sp/s')
	set(gca, 'FontSize', 14)

end
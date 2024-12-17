%% Population Analysis
% J. Fritzinger, updated 1/9/23
clear

% Load in spreadsheet
[base, datapath, ~, ppi] = getPaths();
addpath '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/scripts/helper-functions'
sheetpath = 'scripts/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);


%% Only get sessions with synthetic timbre 

nat(:,1) = cellfun(@(s) contains(s, 'R'), sessions.Bassoon);
nat(:,2) = cellfun(@(s) contains(s, 'R'), sessions.Oboe);
nat(:,3) = cellfun(@(s) contains(s, 'R'), sessions.Other);

any_synth = any(nat(:,1:3), 2);
table = sessions(any_synth, :);

%% Calculate predictable variance for all neurons 

num_index = size(table,1);
V_p = NaN(num_index, 3);
for isesh = 1:num_index

	% Load in session
	putative = table.Putative_Units{isesh};
	CF = table.CF(isesh);
	load(fullfile(datapath, [putative '.mat']))

	% Analysis
	for ispl = 1:3
		params_NT = data(12+ispl, 2);
		if ~isempty(params_NT{1})
			data_ST = analyzeNT(params_NT{1});
			V_p(isesh, ispl) = data_ST.V_p;
		end
	end
end

%% Plot
figure('position', [519,533,1150,305])
tiledlayout(1, 3, 'TileSpacing','tight', 'Padding','compact')
spls = {'Bassoon', 'Oboe', 'Other'};
for ispl = 1:3
	
	nexttile(ispl)
	histogram(V_p(:,ispl), 20)
	%scatter(CFs/1000, V_p(:,ispl), 'filled')
	number = V_p(:,ispl);
	number(isnan(number)) = [];
	title([spls{ispl} ', n=' num2str(length(number))])
	if ispl == 1
		ylabel('# neurons')
	end
	set(gca, 'fontsize', 17)
	grid on
	xlabel('Predictable Variance')
	xlim([0 1])
	ylim([0 100])
end

%% Save figure

savepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/figures/manuscript';
exportgraphics(gcf, fullfile(savepath, 'predictable_variance_NT.png'), 'Resolution', 600)

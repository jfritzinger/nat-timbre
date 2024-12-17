%% plot_temporal_results.m
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
for j = 5

	% Load in data
	putative = table.Putative_Units{j};
	CF = table.CF(j);
	MTF_shape = table.MTF{j};
	load(fullfile(datapath, [putative '.mat']))

	% Analyze data
	param = data(13:14, 2);
	data_NT = cell(2, 1);
	for ii = 1:2
		data_NT{ii} = analyzeNT(param{ii});
	end

	% Plot
	for itype = 1:2
		figure('Position',[865,255,413,692])
		axi = 0;
		num_plots = length(data_NT{itype}.pitch_num);
		for ii = 1:num_plots

			% Set up plots
			ax = axes('Units','normalized',...
				'Position',[0.12 0.05+0.9*axi/num_plots 0.85 0.9/num_plots],...
				'FontSize',6); % [left bottom width height]
			axi = axi + 1;

			% Get data
			j = data_NT{itype}.j(ii,:);
			x = data_NT{itype}.t_spike_rel(ismember(data_NT{itype}.rel_id,j));
			y = data_NT{itype}.rep(data_NT{itype}.rel_id(ismember(data_NT{itype}.rel_id,j)));

			% Plot Histograms
			x = x(x<0.3e6);
			y = y(x<0.3e6);
			plot(x/1000,y,'.k','MarkerSize',10)
			if ii == num_plots
				title([param{itype}.target ' Rasters'])
			end

			% Figure parameters
			if axi == 1
				xlabel('Spike time (ms)');
			else
				ax.XTickLabel = '';
			end
			text(0,0,[num2str(data_NT{itype}.pitch_num(ii)) 'Hz  '],...
				'HorizontalAlignment','right',...
				'VerticalAlignment','bottom', 'FontSize',12)
			ax.FontSize = 16;
			ax.TickLength = [0 0];
			ax.YTickLabel = [];
			ylim([-1 21])
		end
	end


	% for ii = 1:2
	% 	cluster = param{ii}.cluster;
	% 	params = param(ii);
	% 	stims = param{ii}.stims;
	% 	type = 'Rasters';
	% 	[params, fig, ~] = plotPhysNT_Rasters(cluster, params, [], type, data_NT{1});
	% end

end




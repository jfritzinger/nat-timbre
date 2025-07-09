%% supp2_example_neurons

%% Load in 4 examples


target = 'Bassoon';
[base, datapath, ~, ppi] = getPathsNT();
load(fullfile(base, 'model_comparisons', ['Neuron_Time_F0_' target '.mat']), ...
	"neuron_time_F0")

spreadsheet_name = 'Data_Table.xlsx';
sessions = readtable(fullfile(base, spreadsheet_name), ...
	'PreserveVariableNames',true);

accuracy = [neuron_time_F0.accuracy];
[acc_sorted, index] = sort(accuracy, 'descend');

%% Set up figure

figure('Position',[50 50 6.5*ppi, 9*ppi])
%tiledlayout(3, 4, 'TileIndexing','columnmajor')
scattersize = 5;
titlesize = 9;
fontsize = 8;
labelsize = 12;
legsize = 7;
linewidth = 1;
spont_color = [0.4 0.4 0.4];
CF_color = [0.7 0.7 0.7];
rastercolors = {[31,88,239]/256, [218,14,114]/256, [255,176,0]/256};

%% Plots

positions = [1, 2, 3;4, 5, 6;7, 8, 9;10, 11, 12];
for ineuron = 4 %1:4

	% Get best unit
	putative = neuron_time_F0(index(ineuron)).putative;
	load(fullfile(datapath, [putative '.mat']), 'data');

	% Find example in spreadsheet
	s_ind = strcmp(sessions.Putative_Units, putative);
	CF = sessions.CF(s_ind);
	MTF = sessions.MTF(s_ind);
	F0s = getF0s(target);

	% Plot dot raster
	params_NT = data(7, 2);
	data_NT = analyzeNT(params_NT{1});
	temporal = analyzeNT_Temporal(data_NT, CF);

	% Plot dot rasters
	h(positions(ineuron, 1)) = subplot(2, 6, positions(ineuron, 1));
	hold on
	num_stim = 40;
	for j = 1:num_stim
		offset = (j-1)*20; % Adjust offset amount
		if temporal.VS_p(j)<0.01
			scatter(temporal.x{j}/1000,temporal.y{j}+offset,1, 'filled', 'MarkerFaceColor',rastercolors{mod(j, 3)+1})
		else
			scatter(temporal.x{j}/1000,temporal.y{j}+offset,1, 'filled', 'MarkerFaceColor','k', 'MarkerFaceAlpha',0.5)
		end
		yline(offset, 'k')
	end
	ylim([0 20*num_stim])
	xlabel('Time (ms)')
	yticks(linspace(20/2, 20*num_stim-20/2, num_stim))
	yticklabels(round(F0s))
	title('Rasters')
	set(gca,'fontsize',fontsize)
	ylabel('F0 (Hz)')
	xlim([0 0.15])

	% Plot period PSTH
	h(positions(ineuron, 2)) = subplot(2, 6, positions(ineuron, 2));
	hold on
	max_rate = max([temporal.p_hist{:}])-25; %max(temporal.p_hist, [], 'all');
	for j = 1:num_stim

		% Plot PSTHs
		counts = temporal.p_hist{j,:};
		edges = temporal.t_hist{j,:};
		t_bin = edges(1:end-1) + diff(edges)/2; % Bin centers
		x_patch = repelem(edges, 2);
		y_patch = repelem([0; counts(:); 0]', 2);
		y_patch = y_patch(2:end-1); % Creates [0 y1 y1 0 0 y2 y2 0...]
		offset = (j-1)*max_rate; % Adjust offset amount
		if temporal.VS_p(j)<0.01
			patch(x_patch, y_patch + offset, rastercolors{mod(j, 3)+1}, 'FaceAlpha',0.8, 'EdgeColor','k');
		else
			patch(x_patch, y_patch + offset, 'k', 'FaceAlpha',0.3, 'EdgeColor','k');
		end
		period = 1/data_NT.F0s_actual(j)*1000;
		yline(offset, 'k')
		plot([period period], [offset j*max_rate], 'k');
	end
	VS_p = temporal.VS_p;
	r_splithalf(:) = temporal.r_splithalf; % Caclulate reliability metric
	ylim([0 max_rate*num_stim])
	xlabel('Time (ms)')
	xticks(0:30)
	xticklabels({'0', '', '', '', '', '5', '','','','','10','','','','','15', '',''})
	yticks([])
	%yticks(linspace(max_rate/2, max_rate*num_stim-max_rate/2, num_stim))
	xlim([0 17.5])
	yticklabels([])
	xtickangle(0)
	grid on
	title('Period PSTH')
	set(gca,'fontsize',fontsize)


	% Plot VS to harmonics
	h(positions(ineuron, 3)) = subplot(2, 6, positions(ineuron, 3));
	[F0s, peak_harm, peak_harm_num] = calcManualF0s(target);

	VS_harms2 =flipud(temporal.VS_harms);
	peak_harm2 = fliplr(peak_harm_num);
	imagesc(1:30, 1:40, VS_harms2)
	hold on
	for j = 1:num_stim
		rectangle('position', [peak_harm2(j)-0.5 j-0.5, 1, 1], ...
			'EdgeColor','k', 'LineWidth',1)
	end
	xlim([0.51 12.5])
	xlabel('Harmonic #')
	yticklabels([])
	if ineuron == 2 || ineuron == 4
		c = colorbar;
		c.Label.String = 'Vector Strength';
	end
	title('        VS to Harms')
	set(gca,'fontsize',fontsize)

end

%% Arrange figure

left = [linspace(0.06, 0.32, 3) linspace(0.06, 0.32, 3)+0.47];
bottom = [0.04 0.54];
width = 0.12;
height = 0.43;

set(h(1), 'position', [left(1) bottom(2) width height])
set(h(2), 'position', [left(2) bottom(2) width height])
set(h(3), 'position', [left(3) bottom(2) width height])

set(h(4), 'position', [left(4) bottom(2) width height])
set(h(5), 'position', [left(5) bottom(2) width height])
set(h(6), 'position', [left(6) bottom(2) width height])

set(h(7), 'position', [left(1) bottom(1) width height])
set(h(8), 'position', [left(2) bottom(1) width height])
set(h(9), 'position', [left(3) bottom(1) width height])

set(h(10), 'position', [left(4) bottom(1) width height])
set(h(11), 'position', [left(5) bottom(1) width height])
set(h(12), 'position', [left(6) bottom(1) width height])


%% Annotations

labelleft= [0 0.47];
labelbottom = [0.44 0.94];

annotation('textbox',[labelleft(1) labelbottom(2) 0.071 0.058],...
	'String','A','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(2) labelbottom(2) 0.071 0.058],...
	'String','B','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(1) labelbottom(1) 0.071 0.058],...
	'String','C','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(2) labelbottom(1) 0.071 0.058],...
	'String','D','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');


%% Save figure

save_fig = 1;
if save_fig == 1
	filename = 'supp2_pitch_single_unit_example';
	save_figure(filename)
end
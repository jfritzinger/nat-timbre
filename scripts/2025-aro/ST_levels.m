%% ST_level_examples

% Fig___ How responses change/don't change over level
clear

%% Load in spreadsheet

[base, datapath, savepath, ppi] = getPaths();
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(datapath, "data-cleaning/", spreadsheet_name), 'PreserveVariableNames',true);

%% Set up figure
figure('Position',[77,396,695,340])
%tiledlayout(1, 4, 'TileSpacing','compact')
linewidth = 1.5;
fontsize = 20;

%% Plot
for ineuron = 2
	switch ineuron
		case 1
			putative = 'R24_TT2_P13_N03'; % BS
			%putative = 'R24_TT1_P12_N01'; % Low CF
		case 2
			putative = 'R27_TT2_P8_N02'; % BS
			%putative = 'R27_TT3_P1_N01'; % Low CF
		case 3
			putative = 'R24_TT2_P13_N05'; % BS
			%putative = 'R27_TT2_P8_N02'; % Med CF
		case 4
			putative = 'R27_TT3_P7_N08'; % BE
			%putative = 'R27_TT2_P8_N03'; % Med CF
		case 5
			putative = 'R27_TT3_P7_N14'; % BE
			%putative = 'R27_TT4_P8_N10'; % High CF
		case 6
			putative = 'R29_TT3_P5_N10'; % BE
			%putative = 'R25_TT4_P8_N05'; % High CF
	end

	% Load in data
	[base, datapath, savepath, ppi] = getPaths();
	filename = sprintf('%s.mat', putative);
	load(fullfile(datapath,'neural_data', filename)), 'data';
	index = find(cellfun(@(s) strcmp(putative, s), sessions.Putative_Units));
	CF = sessions.CF(index);

	params = data(6:9, 2);
	params = params(~cellfun(@isempty, params));
	data_ST  = analyzeST(params, CF);

	% RM to get spont
	params_RM = data{2, 2};
	data_RM = analyzeRM(params_RM);
	spont = data_RM.spont;

	num_ds = size(data_ST, 2);
	if num_ds == 4
		%data_colors = {'#1b9e77', '#d95f02', '#7570b3', '#e7298a'};
		%data_colors = {'#ffffcc', '#a1dab4', '#41b6c4', '#225ea8'};
		data_colors = {'#82BB95', '#3F985C', '#03882F', '#034E1C'};
	else
		%data_colors = {'#1b9e77', '#d95f02', '#e7298a'};
		%data_colors = {'#a1dab4', '#41b6c4', '#225ea8'};
		data_colors = {'#034E1C', '#03882F', '#82BB95'};
	end

	% Sort
	spls = cell2mat(cellfun(@(p) p.spl, data_ST, 'UniformOutput',false));
	[~, order] = sort(spls);
	order = fliplr(order);
	max_rate = max(cellfun(@(d) max(d.rate), data_ST));
	%nexttile([1 2])
	h(1) = subplot(1, 3, 1);
	hold on
	label_ind = 1;
	for ind = 1:num_ds

		rate = data_ST{order(ind)}.rate;
		rate_std = data_ST{order(ind)}.rate_std;
		rlb = data_ST{order(ind)}.rlb;
		rub = data_ST{order(ind)}.rub;
		fpeaks = data_ST{order(ind)}.fpeaks;
		spl = data_ST{order(ind)}.spl;
		rate_sm = data_ST{order(ind)}.rates_sm;

		% Plot
		rates_sm = smooth_rates(rate, rlb, rub, CF);
		errorbar(fpeaks./1000, rate, rate_std/sqrt(30), 'linestyle', 'none', 'linewidth', 0.8, 'color', data_colors{ind})
		plot(fpeaks./1000, rate, 'LineWidth',linewidth, 'Color',data_colors{:,ind})
		%plot(fpeaks, rate_sm, 'linewidth', linewidth, 'color', data_colors{ind})
		label{1} = [num2str(params{order(ind)}.spl) ' dB SPL'];
		%label_ind = label_ind+1;

		plot_range = [params{1}.fpeaks(1) params{1}.fpeaks(end)];
		xline(CF./1000, '--', 'Color', [0.4 0.4 0.4], 'linewidth', linewidth); % Add CF line
		%label{2} = 'Estimated CF';
		yline(spont, 'k', LineWidth=linewidth)
		%yline(data_RM.spont, 'Color','k', 'LineWidth',2)
		%label(label_ind+1) = {'Spont'};
	end
	plot_range = [params{1}.fpeaks(1) params{1}.fpeaks(end)];
	xline(CF, '--', 'Color', [0.4 0.4 0.4], 'linewidth', linewidth); % Add CF line
	label{label_ind} = 'Estimated CF';
	%yline(data_RM.spont, 'Color','k', 'LineWidth',2)
	%label(label_ind+1) = {'Spont'};
	xlabel('Spectral Peak  Freq. (kHz)')
	ylabel('Avg. rate (sp/s)')
	set(gca, 'Fontsize', fontsize, 'XTick', (plot_range(1)+200:400:plot_range(2)-200)./1000);
	xlim(plot_range./1000);
	grid on
	title('Changes over Level', 'FontSize',fontsize);

end

%%

[base, datapath, savepath, ppi] = getPaths();
tables = readtable(fullfile(datapath,"LMM", "peak_picking_excludeflat.xlsx"));

% Find sessions for target synthetic timbre response
all_neurons = tables.Putative;
neurons = unique(all_neurons);
num_units = size(neurons, 1);
isbin = tables.binmode == 2;
is200 = tables.F0 == 200;

SPLs = [43, 63, 73, 83];
qs = NaN(num_units, 4);
for isesh = 1:num_units

	putative = neurons{isesh};
	isput = cellfun(@(s) strcmp(s, putative), tables.Putative);

	for ispl = 1:4
		ind = isput & isbin & is200 & tables.SPL==SPLs(ispl);
		if any(ind)
			qs(isesh, ispl) = tables.Q(ind);
			qs_log(isesh, ispl) = tables.Q_log(ind);
			CF_group(isesh) = tables.CF_Group(ind);
		end
	end
end


ilow = cellfun(@(d) strcmp('Low', d), CF_group);
imed = cellfun(@(d) strcmp('Med', d), CF_group);
ihigh = cellfun(@(d) strcmp('High', d), CF_group);

spls = [43, 63, 73, 83];
differences = NaN(length(qs), 1);
for ii = 1:length(qs)
	if ~isnan(qs(ii,3)) && ~isnan(qs(ii,1))
		differences(ii) = qs((ii),3)-qs((ii),1);
	elseif ~isnan(qs(ii,2)) && ~isnan(qs(ii,1))
		differences(ii) = qs((ii),2)-qs((ii),1);
	elseif ~isnan(qs(ii,3)) && ~isnan(qs(ii,2))
		differences(ii) = qs((ii),3)-qs((ii),2);
	end
end
decrease = sign(differences)==-1;
increase = sign(differences)==1;

inc = qs(increase,:);
dec = qs(decrease,:);

%nexttile
h(2) = subplot(1, 3, 2);
hold on
for i = 1:4
	swarmchart(ones(length(inc),i)*i, inc(:,i), 15, 'filled')
end
boxplot(inc)
ylim([0 15])
ylabel('Q-value')
xticklabels([43, 63, 73, 83])
xlabel('                            Sound Level (dB SPL)')
title('Sharpening')
set(gca, 'fontsize', fontsize)
plot(1:4, median(inc, 'omitnan'), '--k', 'LineWidth',1.5)

%nexttile
h(3) = subplot(1, 3, 3);
hold on
for i = 1:4
	swarmchart(ones(length(dec),i)*i, dec(:,i), 15, 'filled')
end
boxplot(dec)
ylim([0 15])
yticklabels([])
xticklabels([43, 63, 73, 83])
%xlabel('Sound Level (dB SPL)')
title('Broadening')
set(gca, 'fontsize', fontsize)
plot(1:4, median(dec, 'omitnan'), '--k', 'LineWidth',1.5)


%% Rearrange 

set(h(1), 'position', [0.08 0.17 0.4 0.75])
set(h(2), 'position', [0.55 0.17 0.2 0.75])
set(h(3), 'position', [0.77 0.17 0.2 0.75])

%% Export 

savepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/figures/2025-aro';
set(gcf, 'Renderer', 'painters')
print('-depsc', '-vector', fullfile(savepath,'ST_levels.eps'))



%% Export figure

savepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/figures/2025-aro';
set(gcf, 'Renderer', 'painters')
print('-dsvg', '-vector', fullfile(savepath,'ST_levels.svg'))


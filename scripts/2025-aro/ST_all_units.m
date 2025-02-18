%% ST_all_units 

%% plot_population_all

clear

% Load in spreadsheet
[base, datapath, savepath, ppi] = getPaths();
sheetpath = 'data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(datapath, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);

% Load in spreadsheet with peak information
spreadsheet_name = 'peak_picking.xlsx';
table = readtable(fullfile(datapath, spreadsheet_name));

%% Plot imagesc of all BS responses sorted by CF

spl = [43, 63, 73, 83];
spls = {'43', '63', '73', '83'};
ispl = 2;
MTF_types = {'BS', 'BE', 'H', 'F'};
types = {'Peak', 'Dip', 'Flat'};

for iMTF = 1:3 % iMTF = 1:4

	% Find peaks and dips from table
	isspl = table.SPL == spl(ispl);
	ispeak = strcmp(table.Type, types{iMTF});
	is200 = table.F0 == 200;
	isbin = table.binmode == 2;
	isall = isspl &  ispeak & is200 & isbin;
	putatives = table.Putative(isall);
	peak_freqs = table.Freq(isall);
	num_index = size(putatives, 1);

	CFs = table.CF(isall);
	CF_names = cell(num_index, 1);
	if iMTF == 1
		array_z1 = NaN(num_index,10000);
	elseif iMTF == 2
		array_z2 = NaN(num_index,10000);
	elseif iMTF == 3
		array_z3 = NaN(num_index,10000);
	else
		array_z4 = NaN(num_index,10000);
	end

	for isesh = 1:num_index

		% Load in session
		putative = putatives{isesh};
		CF = CFs(isesh);
		peak_freq = peak_freqs(isesh);
		load(fullfile(datapath, 'neural_data', [putative '.mat']))
		params_ST = data(5+ispl, 2);
		CF_names{isesh} = [num2str(round(CFs(isesh))) ' Hz'];

		% Analysis
		data_ST = analyzeST(params_ST, CF);
		data_ST = data_ST{1};
		params_RM = data{2,2};
		data_RM = analyzeRM(params_RM);
		spont = data_RM.spont;

		% General analysis
		rate = data_ST.rates_sm;
		rate = rate - spont;
		fpeaks = data_ST.fpeaks;
		
		%fpeaks_re_CF = log2(fpeaks/CF);
		if iMTF == 1 || iMTF == 2
			fpeaks_re_CF = log2(fpeaks/peak_freq);
		else
			fpeaks_re_CF = log2(fpeaks/CF);
		end

		% Align by CF (approximately)
		f = linspace(-3, 3, 10000);
		[~, f_ind(1)] = min(abs(fpeaks_re_CF(2)-f));
		[~, f_ind(2)] = min(abs(fpeaks_re_CF(end)-f)); % find indices
		f_interp = linspace(f(f_ind(1)),f(f_ind(2)), f_ind(2)-f_ind(1));

		% Interpolate & get z-score
		r_interp = interp1(fpeaks_re_CF, rate,f_interp, 'spline');
		z_rate = zscore(r_interp);
		%z_rate = r_interp;

		if iMTF == 1
			array_z1(isesh, f_ind(1):f_ind(2)-1) = z_rate;
			CFs1 = CFs;
		elseif iMTF == 2
			array_z2(isesh, f_ind(1):f_ind(2)-1) = z_rate;
			CFs2 = CFs;
		elseif iMTF == 3
			array_z3(isesh, f_ind(1):f_ind(2)-1) = z_rate;
			CFs3 = CFs;
		else
			array_z4(isesh, f_ind(1):f_ind(2)-1) = z_rate;
			CFs4 = CFs;
		end
	end
end


%% Plot

% Set up figure
figure('position', [60,30,950,730])
backgroundcolor =  [0.9 0.9 0.9];
tiledlayout(5, 3, 'TileIndexing','columnmajor')
fontsize = 18;
%titlesize = 18;

for ii = 1:3

	% Order by CF
	if ii == 1
		[~, max_ind] = sort(CFs1);
		CF_order = array_z1(max_ind,:);
		%CFs_ordered = CF_names(max_ind);
	elseif ii == 2
		[~, max_ind] = sort(CFs2);
		CF_order = array_z2(max_ind,:);
	else
		[~, max_ind] = sort(CFs3);
		CF_order = array_z3(max_ind,:);
	end

	% Plot images 
	if ii == 1
		h(1) = subplot(5, 3, [1 4 7 10]);
		%nexttile([5, 1])
	elseif ii == 2
		h(3) = subplot(5, 3, 2);
	else
		h(5) = subplot(5, 3, 3);
		%nexttile
	end
	hh = imagesc(f, 1:size(CF_order, 1), CF_order);
	xline(0, 'k')
	set(hh, 'AlphaData', ~isnan(CF_order))
	set(gca,'color',backgroundcolor);
	if ispl == 1
		ylabel('Neuron Number', 'Color','w')
	end

	yticklabels([])
	xlim([-1 1])
	xticks([-1 0 1])
	clim([-2.2 2.7])
	xticklabels([])

	% Plot overlayed responses
	if ii == 1
		h(2) = subplot(5, 3, 13);
	elseif ii == 2
		h(4) = subplot(5, 3, 5);
	else
		h(6) = subplot(5, 3, 6);
	end
	%nexttile
	hold on
	for iii = 1:size(CF_order, 1)
		patch([f,NaN],[CF_order(iii,:),NaN],'w','EdgeColor','k','LineWidth',1.5,'EdgeAlpha',0.2);
	end
	if ii == 1
		xlabel('Spec. Peak Freq w.r.t. Peak (Oct.)')
	elseif ii == 2
		xlabel('Spec. Peak Freq w.r.t. Dip (Oct.)')
	else
		xlabel('Spec. Peak Freq w.r.t. CF (Oct.)')
	end
	xline(0, 'k')
	yline(0, 'k')
	xlim([-1 1])
	ylim([-2.2 2.7])
	ylabel('Z-score')
	set(gca, 'fontsize', fontsize)
end

% Rearrange
left = linspace(0.05, 0.7, 3);
bottom = linspace(0.07, 0.8, 5);
width = 0.25;
height = 0.155;

set(h(1), 'position', [left(1) 0.23 width 0.73])
set(h(2), 'position', [left(1) bottom(1) width height])

set(h(3), 'position', [left(2) 0.84 width 0.12])
set(h(4), 'position', [left(2) 0.67 width height])

set(h(5), 'position', [left(3) 0.84 width 0.12])
set(h(6), 'position', [left(3) 0.67 width height])

%% Export 

savepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/figures/2025-aro';
set(gcf, 'Renderer', 'painters')
print('-dsvg', '-vector', fullfile(savepath,'ST_all_units1.svg'))


%% Histogram 

figure('Position',[560,667,430,220])

% Plot histogram  
spl = [43, 63, 73, 83];
types = {'Flat', 'Peak', 'Dip'};
isBin = table.binmode == 2;
for ispl = 2
	isSPL = table.SPL == spl(ispl);

	for iMTF = 1:4

		if iMTF == 1
			MTF_target = 'BS';
			isMTF = strcmp(table.MTF, MTF_target);
		elseif iMTF == 2
			MTF_target = 'BE';
			isMTF = strcmp(table.MTF, MTF_target);
		elseif iMTF == 3
			MTF_target = 'Hybrid';
			isMTF = contains(table.MTF, 'H');
		else
			MTF_target = 'F';
			isMTF = strcmp(table.MTF, MTF_target);
		end
		index = isSPL & isMTF & isBin;

		num_dip = sum(cellfun(@(s) strcmp(s, 'Dip'), table.Type(index)));
		num_peak = sum(cellfun(@(s) strcmp(s, 'Peak'), table.Type(index)));
		num_flat = sum(cellfun(@(s) strcmp(s, 'Flat'), table.Type(index)));
		all = sum([num_peak num_dip num_flat]);

		percent_peak(iMTF) = num_peak/all*100;
		percent_dip(iMTF) = num_dip/all*100;
		percent_flat(iMTF) = num_flat/all*100;
	end
	percent_all = [percent_peak; percent_dip; percent_flat]';

	bar(percent_all, 'stacked')
	%title([num2str(spl(ispl)) ' dB SPL'])
	xticklabels({'BS', 'BE', 'Hybrid', 'Flat'})
	legend('Peak', 'Dip', 'Sloping', 'Location','northeastoutside')
	ylabel('% Neurons')
	xlabel('MTF Type')
	ylim([0 100])
	set(gca, 'fontsize', fontsize)

end

%% Export 

savepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/figures/2025-aro';
set(gcf, 'Renderer', 'painters')
print('-dsvg', '-vector', fullfile(savepath,'ST_all_units2.svg'))

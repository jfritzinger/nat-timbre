%%
clear
save_fig = 1;

%% Load in data

[base, ~, ~, ppi] = getPathsNT();
filepath = fullfile(base, 'model_comparisons');
load(fullfile(filepath,  'Model_NT3.mat'), 'nat_model')

% Timbre Rate
load(fullfile(filepath, 'Model_N_Rate_Timbre_All.mat'), "neuron_rate_timbre")
model_rate_timbre = neuron_rate_timbre;
load(fullfile(filepath, 'Neuron_Rate_Timbre_All.mat'), "neuron_rate_timbre")
[~, index] = ismember({model_rate_timbre.putative}, {neuron_rate_timbre.putative});
neuron_rate_timbre = neuron_rate_timbre(index);

% Bassoon Rate
load(fullfile(filepath, 'Model_N_Rate_F0_Bassoon.mat'), "neuron_rate_F0")
model_rate_bassoon = neuron_rate_F0;
load(fullfile(filepath, 'Neuron_Rate_F0_Bassoon.mat'), "neuron_rate_F0")
[~, index] = ismember({model_rate_bassoon.putative}, {neuron_rate_F0.putative});
neuron_rate_bassoon = neuron_rate_F0(index);

% Oboe Rate
load(fullfile(filepath, 'Model_N_Rate_F0_Oboe.mat'), "neuron_rate_F0")
model_rate_oboe = neuron_rate_F0;
load(fullfile(filepath, 'Neuron_Rate_F0_Oboe.mat'), "neuron_rate_F0")
[~, index] = ismember({model_rate_oboe.putative}, {neuron_rate_F0.putative});
neuron_rate_oboe = neuron_rate_F0(index);

% Timbre Timing
load(fullfile(filepath, 'Model_N_Time_Timbre_All.mat'), "neuron_time_timbre")
model_time_timbre = neuron_time_timbre;
load(fullfile(filepath, 'Neuron_Time_Timbre_All.mat'), "neuron_time_timbre")
[~, index] = ismember({model_time_timbre.putative}, {neuron_time_timbre.putative});
neuron_time_timbre = neuron_time_timbre(index);

% Bassoon Timing
load(fullfile(filepath, 'Model_Neuron_Time_F0_Bassoon.mat'), "neuron_time_F0")
model_time_bassoon = neuron_time_F0;
load(fullfile(filepath, 'Neuron_Time_F0_Bassoon.mat'), "neuron_time_F0")
[~, index] = ismember({model_time_bassoon.putative}, {neuron_time_F0.putative});
neuron_time_bassoon = neuron_time_F0(index);

% Oboe Timing
load(fullfile(filepath, 'Model_Neuron_Time_F0_Oboe.mat'), "neuron_time_F0")
model_time_oboe = neuron_time_F0;
load(fullfile(filepath, 'Neuron_Time_F0_Oboe.mat'), "neuron_time_F0")

% Get overlapped indices and cut data down to model size
putative_model = {model_time_oboe.putative};
putative_neuron = {neuron_time_F0.putative};
index = zeros(length(putative_model), 1);
for ind = 1:length(putative_model)
	putative = putative_model{ind};
	possible_ind = find(cellfun(@(p) strcmp(p, putative), putative_neuron));
	if ~isempty(possible_ind)
		index(ind) = possible_ind(1);
	else
		index(ind) = ind;
	end
end
neuron_time_oboe = neuron_time_F0(index);


%% Set up figure

% 4 rows, 7 columns
figure('Position',[50,50,7*ppi,6*ppi])
fontsize = 7;
scattersize = 10;
linewidth = 1;
legsize = 6;

% A.

% Load in spreadsheet
datapath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data';
spreadsheet_name = 'model_r2_values_NT2.xlsx';
sessions = readtable(fullfile(datapath, spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);
isOboe = strcmp(sessions.Instrument, 'Oboe');
isBass = strcmp(sessions.Instrument, 'Bassoon');
edges = linspace(0, 1, 40);


for ii = 1:3
	h(ii) = subplot(5, 7, ii);
	hold on
	edges = linspace(0, 1, 40);
	if ii == 1
		histogram(sessions.SFIE_R2(isOboe), edges, 'FaceAlpha',0.5)
		histogram(sessions.SFIE_R2(isBass), edges, 'FaceAlpha',0.5)
		title('SFIE')
		ylabel('# Neurons')
		ylim([0 120])
		yticks(0:20:120)
	elseif ii == 2
		histogram(sessions.Energy_R2(isOboe), edges, 'FaceAlpha',0.5)
		histogram(sessions.Energy_R2(isBass), edges, 'FaceAlpha',0.5)
		title('Energy')
		xlabel('Variance Explained (R^2)')
		ylim([0 120])
		yticks(0:20:120)
		yticklabels([])
	else
		histogram(sessions.Lat_Inh_R2(isOboe), edges, 'FaceAlpha',0.5)
		histogram(sessions.Lat_Inh_R2(isBass), edges, 'FaceAlpha',0.5)
		title('Broad Inhibition')
		ylim([0 120])
		yticks(0:20:120)
		yticklabels([])
		hleg = legend('Bassoon', 'Oboe');
		hleg.ItemTokenSize = [8,8];
	end
	set(gca, 'Fontsize', fontsize)
	grid on
	xlim([0 0.4])
end


% B.

% Plot
positions = [8, 15, 9, 16, 10, 17];
for ii = 1:6
	h(3+ii) = subplot(5, 7, positions(ii));

	if ii == 1
		data = neuron_rate_timbre;
		model = model_rate_timbre;
		msg = 'Instrument Rate';
	elseif ii == 2
		data = neuron_time_timbre;
		model = model_time_timbre;
		msg = 'Instrument Time';
	elseif ii == 3
		data = neuron_rate_bassoon;
		model = model_rate_bassoon;
		msg = 'Bassoon F0 Rate';
	elseif ii == 4
		data = neuron_time_bassoon;
		model = model_time_bassoon;
		msg = 'Bassoon F0 Time';
	elseif ii == 5
		data = neuron_rate_oboe;
		model = model_rate_oboe;
		msg = 'Oboe F0 Rate';
	else
		data = neuron_time_oboe;
		model = model_time_oboe;
		msg = 'Oboe F0 Time';
	end
	accuracy_data = [data.accuracy];
	accuracy_model = [model.accuracy];

	% Get overlapped indices and cut data down to model size
	putative_model = {model.putative};
	index = NaN(1, length(putative_model));
	for ind = 1:length(model)
		putative = putative_model{ind};
		index(ind) = find(cellfun(@(p) strcmp(p, putative), {nat_model.putative}));
	end
	for ind = 1:length(index)
		if isempty(nat_model(index(ind)).RM_R2)
			nat_model(index(ind)).RM_R2 = 0;
		end
		if isempty(nat_model(index(ind)).MTF_R2)
			nat_model(index(ind)).MTF_R2 = 0;
		end
	end
	RM_R2 = [nat_model(index).RM_R2]';
	MTF_R2 = [nat_model(index).MTF_R2]';
	both_R2 = mean([RM_R2, MTF_R2], 2);

	%good_MTF = MTF_R2 > 0.4;
	%good_MTF = RM_R2 > 0.4;
	good_MTF = both_R2 > 0.5;

	hold on
	scatter(accuracy_data(~good_MTF), accuracy_model(~good_MTF),scattersize, 'filled',...
		'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.5, 'MarkerFaceColor',[0.4 0.4 0.4])
	scatter(accuracy_data(good_MTF), accuracy_model(good_MTF),scattersize, 'filled',...
		'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.8, 'MarkerFaceColor','r')
	if ii == 1 || ii == 2
		xlim([0.4 1])
		ylim([0.4 1])
	else
		xlim([0 0.75])
		ylim([0 0.75])
	end
	if ii == 1 || ii == 2
		ylabel('Model Accuracy')
	end
	if ii == 2 || ii == 4 || ii == 6
		xlabel('Data Accuracy')
	end
	title(msg)
	grid on
	plot([0 1], [0 1], 'k')
	set(gca, 'FontSize', fontsize)

	mdl = fitlm(accuracy_data, accuracy_model);
	p = mdl.Coefficients{2,4};
	fprintf('Instrument, Rate: p=%0.4f\n', p)
	mdl = fitlm(accuracy_data(good_MTF), accuracy_model(good_MTF));
	p = mdl.Coefficients{2,4};
	fprintf('Instrument, Rate (Good): p=%0.4f\n', p)
	[~,p,ci,stats] = ttest(accuracy_data, accuracy_model);
	fprintf('Instrument, Rate: Dist, p=%0.4f\n', p)
	disp('----')

	if ii == 4
	good_MTF = find(both_R2 > 0.5);
	scatter(accuracy_data(good_MTF(21)), accuracy_model(good_MTF(21)), ...
		scattersize, 'blue', 'filled')
	scatter(accuracy_data(good_MTF(22)), accuracy_model(good_MTF(22)),...
		scattersize, 'blue', 'filled')
	end

end

% C.

target = 'Oboe';

load(fullfile(base, 'model_comparisons', 'pop_timing_F0_Oboe_subset2.mat'), ...
	"accur_all","C_all", "num_neurons", "mean_acc", "std_acc")
accur_all2 = accur_all;
C_all2 = C_all;
num_neurons2 = num_neurons;
mean_acc2 = mean_acc;
std_acc2 = std_acc;
load(fullfile(base, 'model_comparisons', ['pop_timing_F0_' target '_subset.mat']), ...
	"accur_all","C_all", "num_neurons")
load(fullfile(base, 'model_comparisons', ['Pop_Rate_F0_' target '.mat']),...
	"pop_rate_F0")
F0s = getF0s(target);

nmodels = length(num_neurons) + length(num_neurons2);
mean_acc1 = mean(accur_all,2);
mean_acc = [mean_acc1(1:12); mean_acc2(1:3); mean_acc1(13:end); mean_acc2(4:6)];
std_acc1 = std(accur_all, [], 2);
std_acc = [std_acc1(1:12); std_acc2(1:3); std_acc1(13:end); std_acc2(4:6)];
num_neurons = [num_neurons(1:12) num_neurons2(1:3) num_neurons(13:24) num_neurons2(4:6)];

h(10) = subplot(5, 7, 22);
errorbar(num_neurons(1:nmodels/2), mean_acc(1:nmodels/2), std_acc(1:nmodels/2)/sqrt(10));
hold on
xlabel('# Neurons in Model')
ylabel('Oboe Accuracy')
ylim([0 1])
grid on
set(gca, 'FontSize', fontsize)
box off
xlim([0 100])

C_diag_all = NaN(nmodels, length(F0s));
iii = 1;
iiii = 1;
for ii = 1:nmodels
	if ~ismember(ii, [13, 14, 15, 28, 29, 30])
		C_diag = NaN(10, length(F0s));
		for irep = 1:10
			C_diag(irep,:) = diag(C_all{iiii,irep});
		end
		iiii = iiii +1;
		C_diag_all(ii,:) = mean(C_diag);
	else
		C_diag_all(ii,:) = diag(C_all2{iii});
		iii = iii +1;
	end
end
h(11) = subplot(5, 7, 23);
y = F0s;
x = num_neurons(1:nmodels/2);
pcolor(y, x, C_diag_all(1:nmodels/2,:), 'EdgeColor','none')
set(gca, 'xscale', 'log')
clim([0 20])
ylabel('# Neurons in Model')
xlabel('F0 (Hz)')
title('Data')
xticks([50 100 250 500 1000 1500])
xticklabels([50 100 250 500 1000 1500]/1000)
set(gca, 'FontSize', fontsize)
box off

load(fullfile(base, 'model_comparisons', ['Model_Pop_Timing_F0_' target '.mat']), ...
	"accur_all","C_all", "num_neurons")
load(fullfile(base, 'model_comparisons', ['Model_Pop_Rate_F0_' target '.mat']),...
	"pop_rate_F0")

nmodels = length(num_neurons);
mean_acc = mean(accur_all,2);
std_acc = std(accur_all, [], 2);
for ii = 1:nmodels
	C_diag_all(ii,:) = diag(C_all{ii});
end
errorbar(h(10), num_neurons(1:nmodels/2), mean_acc(1:nmodels/2), std_acc(1:nmodels/2)/sqrt(10));
hleg = legend(h(10), 'Data', 'Model', 'Location','southeast');
hleg.ItemTokenSize = [8,8];


h(12) = subplot(5, 7, 24);
y = F0s;
x = num_neurons(1:nmodels/2);
pcolor(y, x, C_diag_all(1:nmodels/2,:), 'EdgeColor','none')
set(gca, 'xscale', 'log')
clim([0 20])
xlabel('F0 (Hz)')
title('Model')
yticks([])
a=colorbar;
a.Label.String = '# Accurate Predictions';
xticks([50 100 250 500 1000 1500])
xticklabels([50 100 250 500 1000 1500]/1000)
set(gca, 'FontSize', fontsize)
box off


% D.

% Choose neuron
target = 'Bassoon';
F0s = calcManualF0s(target);
positions = [28 5 12 19 26 11 18 25; 31 7 14 21 28 13 20 27];
h_pos = [20 14 15 16; 23 18 19 20];
for ii = 1:2
	if ii == 1
		putative = 'R29_TT4_P2_N04'; % good for both
		%putative = 'R29_TT4_P1_N04';
	else
		putative = 'R29_TT4_P1_N02'; % bad for model
		%putative = 'R27_TT3_P1_N05';
	end
	%putative = 'R29_TT3_P1_N02';
	%putative = 'R29_TT3_P2_N11';
	%putative = 'R29_TT4_P2_N06';
	%putative = 'R29_TT3_P3_N02';
	%putative =  'R27_TT2_P11_N03';
	s_ind = strcmp({nat_model.putative}, putative);
	CF = nat_model(s_ind).CF;

	% Two examples
	bass_psth_all = nat_model(s_ind).bass_PSTH_all;
	t = 0:0.25:300/1000; %linspace(0, 300, 1200);
	t = t(1:end-1);
	x = cell(40, 1);
	y = cell(40, 1);
	for istim = 1:40

		bass_psth_one = bass_psth_all{istim};
		spike_rate = [];
		spike_rep = [];
		for irep = 1:20
			rate_fh = bass_psth_one(irep,:);
			T = 0.3;
			EventTimes = genNHPP(rate_fh,T, 1);
			spike_rate = [spike_rate EventTimes];
			num_spikes = length(EventTimes);
			spike_rep = [spike_rep; irep*ones(num_spikes, 1)];
		end
		x{istim} = spike_rate;
		y{istim} = spike_rep;
		spike_rate = spike_rate*1000;

		% Calculate period histogram
		dur = 300;
		onset = 25; % 25 ms onset
		ind_onset = spike_rate<onset;
		spike_times = spike_rate(~ind_onset);
		reps = spike_rep(~ind_onset);
		period = 1000 / (F0s(istim));	% Period in ms
		num_periods = floor((dur-onset)/period);	% Number of full periods in the stimulus
		ind_full_period = spike_times<(num_periods*period);
		subset_spike_times = spike_times(ind_full_period);
		subset_reps = reps(ind_full_period);
		wrapped_times = mod(subset_spike_times, period);
		num_bins = 30; %round(period/0.25); % 30; % Number of bins for histogram
		edges = linspace(0, period, num_bins+1); % Bin edges
		counts = histcounts(wrapped_times, edges); % Create histogram
		temporal.p_hist{istim,:} = counts;
		temporal.t_hist{istim,:} = edges;

		% Calculate VS to harmonics
		for iharm = 1:30
			harm = (F0s(istim)+(F0s(istim)*(iharm-1)));
			period = 1000 / harm; % Get period of each harmonic
			phases = 2 * pi * (mod(subset_spike_times, period) / period);
			VS_harms(iharm) = abs(mean(exp(1i * phases)));
			if ~isempty(phases)
				p_value_harms(iharm) = circ_rtest(phases); % Rayleigh statistic (P < 0.01)
			else
				p_value_harms(iharm) = NaN;
			end
			harms(iharm) = harm;
		end
		temporal.VS_harms(istim,:) = VS_harms;

	end
	load(fullfile(datapath, 'Neural_Data', [putative '.mat']), 'data');

	% Plot RM
	params_RM = data{2,2}; % Gets binaural response map
	data_RM = analyzeRM(params_RM);
	spont_color = [0.4 0.4 0.4];
	for ind = 1:3
		h(h_pos(ii, 1)+ind) = subplot(5, 7, positions(ii, 1)+ind);
		hold on
		if ind == 1
			area(data_RM.freqs,data_RM.rates(:,3),'EdgeColor', '#A49BD0',...
				'FaceColor', '#A49BD0','LineWidth',linewidth) % 33 dB
			spl = '33 dB SPL';
			plot(nat_model(s_ind).freqs, nat_model(s_ind).RM_rate(:,3), 'linewidth', 1.5)
		elseif ind == 2
			area(data_RM.freqs,data_RM.rates(:,4),'EdgeColor', '#5E50A9',...
				'FaceColor', '#5E50A9','LineWidth',linewidth) % 53 dB
			spl = '53 dB SPL';
			plot(nat_model(s_ind).freqs, nat_model(s_ind).RM_rate(:,4), 'linewidth', 1.5)
		elseif ind == 3
			area(data_RM.freqs,data_RM.rates(:,5),'EdgeColor', '#20116B',...
				'FaceColor', '#20116B','LineWidth',linewidth) % 73 dB
			spl = '73 dB SPL';
			plot(nat_model(s_ind).freqs, nat_model(s_ind).RM_rate(:,5), 'linewidth', 1.5)
		end
		plot(data_RM.freqs([1 end]),[1 1]*data_RM.spont,'-','LineWidth',...
			linewidth, 'Color',spont_color)
		text(0.05, 0.90, spl, 'Units', 'normalized', ...
			'VerticalAlignment', 'top', 'FontSize',legsize)
		grid on
		hold off
		ylim([0 max(data_RM.rates, [], 'all')+5])
		yticks([0 40])
		xticks([200 500 1000 2000 5000 10000])
		xticklabels([200 500 1000 2000 5000 10000]./1000)
		set(gca, 'XScale', 'log','fontsize',fontsize);
		if ind == 1
			xlabel('Tone Freq. (kHz)')
			hLegend = legend('', '', '', 'Location','northeast ', ...
				'Box','off');
			hLegend.ItemTokenSize = [12,6];
		elseif ind == 2
			ylabel('Avg. Rate (sp/s)')
			xticklabels([])
		elseif ind == 3
			title('Response Map')
			xticklabels([])
		end
	end


	% Plot MTF
	h(h_pos(ii, 2)) = subplot(5, 7, positions(ii, 2));
	hold on
	params_MTF = data{3,2}; % Gets binaural MTFN
	data_MTF = analyzeMTF(params_MTF);
	line([1 data_MTF.fms(end)], [1 1]*data_MTF.rate(1),'Color',[0.4 0.4 0.4],...
		'LineWidth',linewidth);
	errorbar(data_MTF.fms,data_MTF.rate, ...
		data_MTF.rate_std/sqrt(params_MTF.nrep),'.k', 'LineWidth',...
		0.8, 'CapSize',4);
	line(data_MTF.fms,data_MTF.rate,'Color','k', 'Marker','.', ...
		'MarkerSize',5, 'MarkerFaceColor','w', 'LineWidth', linewidth);
	zeroed_data = data_MTF.rate-data_MTF.rate(1); % Scale model 
	zeroed_model = nat_model(s_ind).MTF_rate-nat_model(s_ind).MTF_rate(1);
	max_data = max(zeroed_data);
	model_rate = zeroed_model/max(zeroed_model)*max_data+data_MTF.rate(1);
	%plot(nat_model(s_ind).fms(2:end), nat_model(s_ind).MTF_rate(2:end), 'k')
	plot(nat_model(s_ind).fms(2:end), model_rate(2:end), 'r', 'linewidth', 1.5)
	set(gca, 'xscale', 'log')
	xticks([2 5 10 20 50 100 200 500])
	xlim([1.2 500])
	set(gca,'fontsize',fontsize)
	hleg = legend('Unmod.', 'Data', '', 'Model', 'location', 'best', 'box', 'off');
	hleg.ItemTokenSize = [8, 8];
	title('MTF')
	ylim([0 110])
	grid on

	% Plot dot rasters
	h(h_pos(ii, 3)) = subplot(5, 7, positions(ii, 3:5));
	rastercolors = {[31,88,239]/256, [218,14,114]/256, [255,176,0]/256};
	num_stim = 40;
	hold on
	for j = 1:num_stim

		x1 = x{j};
		y1 = y{j};

		offset = (j-1)*20; % Adjust offset amount
		scatter(x1, y1+offset,1, 'filled', 'MarkerFaceColor',rastercolors{mod(j, 3)+1})
		yline(offset, 'k')
	end
	ylim([0 20*num_stim])
	xlabel('Time (ms)')
	yticks(linspace(20/2, 20*num_stim-20/2, num_stim))
	yticklabels([])
	grid on
	title('Model')
	xlim([0 0.15])
	set(gca,'fontsize',fontsize)

	% Plot real
	h(h_pos(ii, 4)) = subplot(5, 7, positions(ii, 6:8));
	hold on
	
	params_NT = data(7, 2);
	data_NT = analyzeNT(params_NT{1});
	temporal = analyzeNT_Temporal(data_NT, CF);
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
	set(gca,'fontsize',fontsize)
	ylabel('Bassoon F0 (Hz)')
	xlim([0 0.15])
	title('Data')
end

% Arrange 

left = [0.05 0.18 0.31 0.49 0.61 0.76 0.88];
bottom = [0.06 0.31 0.53 0.79];
width = 0.1;
height = 0.16;
height2 = 0.64;

set(h(1), 'position', [left(1) bottom(4) width height])
set(h(2), 'position', [left(2) bottom(4) width height])
set(h(3), 'position', [left(3) bottom(4) width height])

set(h(4), 'position', [left(1) bottom(3) width height])
set(h(5), 'position', [left(1) bottom(2) width height])
set(h(6), 'position', [left(2) bottom(3) width height])
set(h(7), 'position', [left(2) bottom(2) width height])
set(h(8), 'position', [left(3) bottom(3) width height])
set(h(9), 'position', [left(3) bottom(2) width height])

set(h(10), 'position', [left(1) bottom(1) width-0.01 height])
set(h(11), 'position', [left(2)+0.02 bottom(1) width-0.02 height])
set(h(12), 'position', [left(3)-0.02 bottom(1) width-0.02 height])

%set(h(13), 'position', [left(4) bottom(4) width height])
set(h(14), 'position', [left(5) bottom(4) width height])
set(h(16), 'position', [left(4) bottom(1) width height2])
set(h(15), 'position', [left(5) bottom(1) width height2])

%set(h(17), 'position', [left(6) bottom(4) width height])
set(h(18), 'position', [left(7) bottom(4) width height])
set(h(20), 'position', [left(6) bottom(1) width height2])
set(h(19), 'position', [left(7) bottom(1) width height2])

set(h(21), 'position', [left(4) bottom(4) width height/3])
set(h(22), 'position', [left(4) bottom(4)+height/3 width height/3])
set(h(23), 'position', [left(4) bottom(4)+height/3*2 width height/3])
set(h(24), 'position', [left(6) bottom(4) width height/3])
set(h(25), 'position', [left(6) bottom(4)+height/3 width height/3])
set(h(26), 'position', [left(6) bottom(4)+height/3*2 width height/3])


% Annotate 

labelleft= left-0.04;
labelbottom = bottom+0.16;
labelsize = 12;

annotation('textbox',[labelleft(1) labelbottom(4) 0.071 0.058],...
	'String','A','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(1) labelbottom(3) 0.071 0.058],...
	'String','B','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(1) labelbottom(1) 0.071 0.058],...
	'String','C','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(4) labelbottom(4) 0.071 0.058],...
	'String','D','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(6) labelbottom(4) 0.071 0.058],...
	'String','E','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');

%% Save figure

if save_fig == 1
	filename = 'fig11_model';
	save_figure(filename)
end

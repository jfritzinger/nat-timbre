%% testing_spike_generator
clear

%% Load in data

[base, datapath, savepath, ppi] = getPathsNT();
load(fullfile(base, 'model_comparisons',  'Model_NT.mat'), 'nat_model')

%% Get single neuron

target = 'Bassoon';
F0s = calcManualF0s(target);

% Choose neuron
%putative = 'R29_TT3_P1_N02';
%putative = 'R29_TT3_P2_N11';
%putative = 'R29_TT4_P1_N04';
%putative = 'R29_TT4_P2_N06';
%putative = 'R29_TT3_P3_N02';

s_ind = strcmp({nat_model.putative}, putative);
CF = nat_model(s_ind).CF;

% Two examples
bass_psth_all = nat_model(s_ind).bass_PSTH_all;
t = 0:0.25:300/1000; %linspace(0, 300, 1200);
t = t(1:end-1);

%   Generates a NHPP, N(t): # of events by time t
% INPUTS
%   rate_fh : function handle for rate (vectorized)
%   T : end of time horizon
%   n : number of sample paths to generate (default n = 1)
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
		
		% % Plot
		% figure
		% plot(t, bass_psth_one(1,:))
		% hold on
		% for ii = 1:length(EventTimes)
		% 	xline(EventTimes(ii)*1000)
		% end
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


%% Plot dot rasters

figure
tiledlayout(1, 3, 'TileSpacing','compact')
fontsize = 8;

nexttile
rastercolors = {[31,88,239]/256, [218,14,114]/256, [255,176,0]/256};
num_stim = 40;
hold on
for j = 1:num_stim

	x1 = x{j};
	y1 = y{j};

	offset = (j-1)*20; % Adjust offset amount
	scatter(x1, y1+offset,3, 'filled', 'MarkerFaceColor',rastercolors{mod(j, 3)+1})
	yline(offset, 'k')
end
ylim([0 20*num_stim])
xlabel('Time (ms)')
yticks(linspace(20/2, 20*num_stim-20/2, num_stim))
yticklabels(round(F0s))
grid on
title(CF)

% Plot period histogram

nexttile
hold on
max_rate = max([temporal.p_hist{:}]); %-25; %max(temporal.p_hist, [], 'all');
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
	patch(x_patch, y_patch + offset, rastercolors{mod(j, 3)+1}, 'FaceAlpha',0.8, 'EdgeColor','k');
	period = 1/F0s(j)*1000;
	yline(offset, 'k')
	plot([period period], [offset j*max_rate], 'k');
end
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

nexttile

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
xlabel('Harmonic Number')
yticklabels([])
c = colorbar;
c.Label.String = 'Vector Strength';
title('        VS to Harmonics')
set(gca,'fontsize',fontsize)
clim([0 0.8])


%% Calc ISI
% for j = 1:num_stim
% 
% 	x1 = x{j}';
% 	y1 = y{j};
% 
% 	ISI = arrayfun(@(ii) diff(x1(y1==ii)), 1:20, 'UniformOutput', false);
% 	ISI_all = vertcat(ISI{:});
% 
% 	% figure
% 	% edges = 0:0.001:0.010;
% 	% histogram(ISI_all, edges)
% 	% xlim([0 0.01])
% 	% hold on
% 	% xline(0.001)
% 
% 	num_spikes = length(ISI_all);
% 	perc_less_1ms = sum(ISI_all < 0.001)/num_spikes*100;
% 	disp(['Percent < 1 ms = ' num2str(round(perc_less_1ms))])
% 
% end

%% Plot PSTH
% for istim = 29:40
% 	bass_psth_one = bass_psth_all{istim};
% 	rate_fh = bass_psth_one(1,:);
% 
% 	spikes = x{istim}; % ms
% 	spikereps = y{istim};
% 	min_dis = 0.25;
% 	edges = (0:min_dis:300)/1000;
% 	t = 0+min_dis/2:min_dis:300-min_dis/2;
% 	for irep = 1:20
% 		h_bass(irep, :) = histcounts(spikes(spikereps==irep), edges);
% 	end
% 	h_all = sum(h_bass, 1);
% 
% 
% 	figure
% 	hold on
% 	plot(t, rate_fh./(max(rate_fh)))
% 	plot(t,h_all./max(h_all))
% 	ylim([0 1])
% end

% Start getting a DC offset at 293 Hz F0
% Still entrainment to the spikes at 554 Hz
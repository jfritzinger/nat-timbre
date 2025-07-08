%% testing_spike_generator

%% Load in data

[base, datapath, savepath, ppi] = getPathsNT();
load(fullfile(base, 'model_comparisons',  'Model_NT.mat'), 'nat_model')

%% Get single neuron

target = 'Bassoon';
F0s = getF0s(target);
F0s = log10(F0s);

% Two examples
bass_psth_all = nat_model(2).bass_PSTH_all;
CF = nat_model(4).CF;
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
end


% Plot dot rasters

rastercolors = {[31,88,239]/256, [218,14,114]/256, [255,176,0]/256};
num_stim = 40;
figure
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

% Calc ISI

for j = 1:num_stim

	x1 = x{j}';
	y1 = y{j};

	ISI = arrayfun(@(ii) diff(x1(y1==ii)), 1:20, 'UniformOutput', false);
	ISI_all = vertcat(ISI{:});

	% figure
	% edges = 0:0.001:0.010;
	% histogram(ISI_all, edges)
	% xlim([0 0.01])
	% hold on
	% xline(0.001)

	num_spikes = length(ISI_all);
	perc_less_1ms = sum(ISI_all < 0.001)/num_spikes*100;
	disp(['Percent < 1 ms = ' num2str(round(perc_less_1ms))])

end

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
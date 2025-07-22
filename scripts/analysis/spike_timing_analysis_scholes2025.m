%% Try spike statistics on natural timbre dataset 
%% ------------------------------------------------------------------------
%                    Calculating VS, mode-locking, CI .... etc
% -------------------------------------------------------------------------
clear;
warning('off','all')

% Add the paths required.
%addpath functions;
%addpath functions/SACcode;  % You need a path to the correlation index code. Code is from Joris.

%% Load in data and arrange in cell arrays 

% Load putative unit
[base, datapath, ~, ppi] = getPathsNT();
sessions = readtable(fullfile(base, 'Data_Table.xlsx'), ...
	'PreserveVariableNames',true);
putative = 'R29_TT2_P3_N03';
load(fullfile(datapath, [putative '.mat']), 'data');
s_ind = strcmp(sessions.Putative_Units, putative);
CF = sessions.CF(s_ind);
MTF = sessions.MTF(s_ind);
params_NT = data(7, 2);
data_NT = analyzeNT(params_NT{1});
temporal = analyzeNT_Temporal(data_NT, CF);

% Arrange in spike times into cell arrays 
spike_times = temporal.spike_times';

% Calculate F0s
F0s = calcManualF0s('Bassoon');

%% Options - which statistics to recalculate.

calcVS = true;     % Calculate vector strength statistics.
calcCI = true;     % Calculate correlation index statistics.
calcML = true;     % Calculate mode-locking statistics.
calcGC = true;     % Calculate reliability and envelope fluctuation statistics from Gai and Carney.
calcSACpeaks = true; % Calculate the SAC peaks from spike trains.

%% -----------------------------------------------------------------------
%         This calculates statistics for every condition, every unit.

win = [25 300];
numPHbins = 52;
phbin = 0:1/numPHbins:1; % see Malone et al. 2007
coinc_win = .05;         % Coincidence window for Correlation Index.
duration = 300;
ndata = length(spike_times);
SACfigdir = 'SACfigs';

% Variables for the peak analysis.
global CRITERION
CRITERION =.05;
SACpeaks_minspkpersweep = 6;    % We arbitrarily don't run the SAC peak picking
% % on anything less than 4 spikes per
% % sweep.

% Initialise outputs for which ever are being calculated.
if calcVS
	VS = cell(size(spike_times));
end
if calcML
	ML = cell(size(spike_times));
end
if calcCI
	CI = cell(size(spike_times));
end
if calcGC
	GC = cell(size(spike_times));
end
if calcSACpeaks
	SACpeaks = cell(size(spike_times));
end

completedvector = cell(size(spike_times));
timerVal = tic;
parfor ii = 1:ndata
	timerVal2 = tic;
	theseSpikeTimes = spike_times{ii};
	fprintf('Data %d /%d:\n',ii,ndata);

	if ~isempty([theseSpikeTimes{:}])

		% Select out the right bit of the spike train.
		spktrainset = cellfun(@(x)  x(x>=win(1) & x<=win(2)),theseSpikeTimes,'UniformOutput',false);

		% Some info about the spike train.
		nreps = 20;
		thisF0 = F0s(ii);
		per = 1000 / thisF0;

		% ------------- vector strength ----------------
		if calcVS
			% fprintf('VS,\n');

			% calculate spike-phase as a fraction of the period.
			st = [spktrainset{:}];
			p_st = (st - floor(st/per)*per)/per; % spike times in periods
			l_st = sum(~isnan(p_st));            % number of spikes without nans.

			VS{ii}.VSvalue = sqrt(sum(cos(2*pi*p_st))^2+sum(sin(2*pi*p_st))^2)/l_st;
			VS{ii}.Rayleigh = 2*l_st*(VS{ii}.VSvalue^2);
			VS{ii}.PH = hist(p_st,phbin);
			VS{ii}.PHbins = phbin;
			VS{ii}.spikes_per_rep = length(st)/nreps;
			completedvector{ii}.value = 1;
		end

		% ----------- mode-locking -------------------
		if calcML
			% fprintf('ML,\n');
			% Call getModelockCriteria to perform the full set of tests.
			ML{ii} =getModelockCriteria(theseSpikeTimes, thisF0, win, 1);
			ML{ii}.modFreq = thisF0;
		end

		% ------------- correlation index ----------------

		if calcCI
			% Calcuate the correlation index.
			% fprintf('CI\n');
			maxlagfactor = 1.5;
			maxlag_ms= maxlagfactor*1e3/F0s(ii);
			[h, bc] = SPTCORR(spktrainset,'nodiag',maxlag_ms,coinc_win,duration,'LouageNorm');
			%[h, bc] = SPTCORR(spktrainset,'nodiag',20,coinc_win,duration,'LouageNorm');
			CI{ii}.unsmoothedSAC = h;
			CI{ii}.CIvalue_0lag = h(bc==0); % CI is nominally the  zero lag value.
			CI{ii}.CIvalue_max = max(h); % But that might not be the maximum value.
			CI{ii}.lags = bc;

			% Calculate significance p-value associated with the CI,
			% using the Joris method.
			if sum(~cellfun('isempty', spktrainset))>0
				% remove any empty cells
				spktrainset = spktrainset(~cellfun('isempty', spktrainset));
				pval = sacpeaksign(spktrainset,coinc_win,500,duration);  % Very slow.
				CI{ii}.p = pval;
			else
				CI{ii}.p =  NaN;
			end
		end

		% figure
		% plot(CI{ii}.lags, CI{ii}.unsmoothedSAC)

		% ------------- SAC peak statistics ----------------

		% An elaborate analysis of the SAC which looks for peaks additional
		% to the expected ones fom phase-locking. VERY slow.
		% Cannot emphasize enough that running this is a big deal.
		if calcSACpeaks
			% Calculate the correlation index.
			fprintf('SAC peaks...\n');
			% [SACpeaks{ii}, details] = getSACpeaks_fast(spktrainset,CI{ii},win, ...
			% 	thisF0, SACpeaks_minspkpersweep,SACfigdir,ii);
			[SACpeaks{ii}, details] = getSACpeaks(spktrainset,CI{ii},win, ...
				thisF0, SACpeaks_minspkpersweep,SACfigdir,ii);
			if isfield(SACpeaks{ii},'amfreq')
				SACpeaks{ii}.details = details;
			end
		end

		% ------------- Gai & Carney statistics ----------------

		% Calculates: "reliability" (a PSTH correlation measure) and
		% "envelope fluctuation".
		if calcGC
			GC{ii} = getGaiCarneyStats(spktrainset,win);
		end

		% figure
		% plot(SACpeaks{ii}.lagaxis, SACpeaks{ii}.smoothsac)
		% hold on
		% if ~isempty(SACpeaks{ii}.setofsigpeaklags_p001)
		% 	xline(SACpeaks{ii}.setofsigpeaklags_p001)
		% end
	end
	fprintf('F0 took %0.2f seconds\n', toc(timerVal2))
end
fprintf('Neuron took %0.2f seconds\n', toc(timerVal))

% Save the calculated statistics.
%save datafiles\SpikeStats_tmp;

%% Plot rasters and period histogram

figure
tiledlayout(1, 3)
rastercolors = {[31,88,239]/256, [218,14,114]/256, [255,176,0]/256};

% Plot dot rasters
nexttile
hold on
num_stim = 40;
for j = 1:num_stim
	offset = (j-1)*20; % Adjust offset amount
	if temporal.VS_p(j)<0.01
		scatter(temporal.x{j}/1000,temporal.y{j}+offset,3, 'filled', 'MarkerFaceColor',rastercolors{mod(j, 3)+1})
	else
		scatter(temporal.x{j}/1000,temporal.y{j}+offset,3, 'filled', 'MarkerFaceColor','k', 'MarkerFaceAlpha',0.9)
	end
	yline(offset, 'k')
end
ylim([0 20*num_stim])
xlabel('Time (ms)')
yticks(linspace(20/2, 20*num_stim-20/2, num_stim))
yticklabels(round(F0s))
grid on
xlim([0 0.15])

% Plot period histogram
nexttile
hold on
max_rate = max([temporal.p_hist{:}])-max([temporal.p_hist{:}])/3; %-150; %max(temporal.p_hist, [], 'all');
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
		patch(x_patch, y_patch + offset, rastercolors{mod(j, 3)+1}, 'FaceAlpha',0.9, 'EdgeColor','k');
	else
		patch(x_patch, y_patch + offset, 'k', 'FaceAlpha',0.9, 'EdgeColor','k');
	end
	period = 1/data_NT.pitch_num(j)*1000;
	yline(offset, 'k')
end
VS_p = temporal.VS_p;
ylim([0 max_rate*num_stim])
xlabel('Time (ms)')
xticks(0:30)
yticks(linspace(max_rate/2, max_rate*num_stim-max_rate/2, num_stim))
yticklabels(round(VS_p, 3))
grid on
title('Period PSTH')

nexttile
maxVal = 8;
for j = 1:num_stim
	offset = maxVal*(j-1); % Adjust offset amount
	plot(SACpeaks{j}.smoothsac+offset, 'k')
	hold on
	if ~isempty(SACpeaks{j}.setofsigpeaks)
		correct_lags = SACpeaks{j}.setofsigpeaks;
		for ii = 1:length(correct_lags)
			plot([correct_lags(ii) correct_lags(ii)], [offset offset+maxVal], 'r')
		end
	end
end
% nexttile
% %maxPerCell = cellfun(@(s) max(s.smoothsac), SACpeaks);
% maxValue = 8; %max(maxPerCell);
% for j = 1:num_stim
%     offset = maxValue*(j-1);
%     plot(SACpeaks{j}.lagaxis, SACpeaks{j}.smoothsac+offset, 'k')
%     hold on
%     if ~isempty(SACpeaks{j}.setofsigpeaklags)
%         correct_lags = SACpeaks{j}.setofsigpeaklags;
%         % Align the marker to the closest lagaxis value for display
%         for ii = 1:length(correct_lags)
%             plot([correct_lags(ii) correct_lags(ii)], [offset offset+maxValue], 'r')
%         end
%     end
%     yline(offset)
% end


% nexttile
% for j = 1:num_stim
% 	offset = (j-1); % Adjust offset amount
% 	plot(SACpeaks{j}.lagaxis, SACpeaks{j}.smoothsac+offset, 'k')
% 	hold on
% 	if ~isempty(SACpeaks{j}.setofsigpeaklags)
% 		correct = min(abs(SACpeaks{j}.setofsigpeaklags));
% 		correct_lags = SACpeaks{j}.setofsigpeaklags+correct;
% 		for ii = 1:length(correct_lags)
% 			plot([correct_lags(ii) correct_lags(ii)], [offset offset+1], 'r')
% 		end
% 	end
% 	% plot([0 0], [offset offset+1],'r')
% 	% plot([1000/F0s(j) 1000/F0s(j)], [offset offset+1],'r')
% 	% plot([-1*1000/F0s(j) -1*1000/F0s(j)], [offset offset+1],'r')
% 	yline(offset)
% end



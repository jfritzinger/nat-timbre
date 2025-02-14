%% plot_rvf_prediction_example
clear 

%% Load in natural timbre file names

if ismac
	fpath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data';
else
	fpath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\';
end
tuning = readtable(fullfile(fpath, 'Tuning.xlsx')); % Load in tuning
target = 'Bassoon';

% Get all .wav files containing the target instrument name
listing = dir(fullfile(fpath, 'waveforms', ['*' target '*.wav']));
files = {listing.name};

% Extract note names and find corresponding frequencies
note_names = extractBetween(files, 'ff.', '.');
[~, index] = ismember(note_names, tuning.Note);
F0s = tuning.Frequency(index);

% Sort files and frequencies by pitch
[F0s, order] = sort(F0s);
files = files(order);

% Initialize variables for later use (if needed)
nfiles = numel(files);
wav_npts = zeros(1, nfiles);
wav_data = cell(1, nfiles);

%% Analyze natural timbre

[base, datapath, ~, ppi] = getPaths();
sheetpath = 'data/2025-manuscript/data-cleaning';
spreadsheet_name = 'PutativeTable2.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);

% Find all sessions with RVF 
%hasRVF = cellfun(@(s) contains(s, 'R'), sessions.RVF);

% Example session
putative = 'R29_TT3_P5_N02';
filename = sprintf('%s.mat', putative);
load(fullfile(datapath,'neural_data', filename)), 'data';
index = find(cellfun(@(s) strcmp(putative, s), sessions.Putative_Units));
CF = sessions.CF(index);
MTF_shape = sessions.MTF{index};

% Plot NT
params_NT = data(14, 2); % 13 is oboe, 14 is bassoon
data_NT = cell(2, 1);
if ~isempty(params_NT{1})
	data_NT = analyzeNT(params_NT{1});
end


%% Plot neurograms

% if ismac
% 	savepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/AN_neurogram';
% else
% 	savepath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\AN_neurogram';
% end
if ismac
	savepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/decomped';
else
	savepath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\decomped';
end


vel_CF = NaN(1, nfiles);
for ii = 1:nfiles

	% Get single stimulus
	target = extractBefore(files{ii}, '.');
	target_F0 = extractBetween(files{ii}, 'ff.', '.wav');
	% savefile = sprintf('%s_F0_%s_Velocity.mat', target, target_F0{1});
	% load(fullfile(savepath, savefile),'vel_all', 'CFs_all')

	savefile = sprintf('%s_F0_%d_Velocity.mat', target, round(F0s(ii)));
	load(fullfile(savepath, savefile), 'v', 'v_harms')
	vel_all = v;
	CFs_all = v_harms;
	

	% Find 
	[~, CF_ind] = min(abs(CF-CFs_all));
	vel_CF(ii) = vel_all(CF_ind);

end

%% Get rate / interpolated rate from RVF

% Plot RVF
params_RVF = data{5, 1};
[fig, data_RVF] = plotPhysRVF([], params_RVF, []);

%% Rearrange
[sorted_vel, ind] = sort(data_RVF.velocities);
sorted_rate = data_RVF.rate(ind);

% Interpolate (linear)
RVF_vel = -9:0.1:9;
RVF_rates = interp1(sorted_vel, sorted_rate, RVF_vel, 'linear');

% Match each velocity to a rate 
for ii = 1:nfiles

	% Find rate most similar to velocity
	this_vel = vel_CF(ii);
	[~, index] = min(abs(this_vel-RVF_vel));
	if index == 1 
		rates_vel(ii) = RVF_rates(index+1);
	elseif index == length(RVF_vel)
		rates_vel(ii) = RVF_rates(index-1);
	else
		rates_vel(ii) = RVF_rates(index);
	end
end

%% Plot both 

figure
tiledlayout(1, 3)

nexttile
plot(F0s, vel_CF)
title('Velocities for each stimulus')
xlabel('Pitch (Hz)')
ylabel('Velocity (kHz/ms)')
set(gca, 'xscale', 'log')
set(gca, 'fontsize', 16)
xticks([60 110 220 440])


nexttile
plot(RVF_vel, RVF_rates);
title('Interpolated RVF')
set(gca, 'fontsize', 16)
xlabel('Velocity (kHz/ms)')
ylabel('Avg. Rate (sp/s)')


% Plot natural timbre
nexttile
hold on
errorbar(F0s,data_NT.rate, data_NT.rate_std/sqrt(params_NT{1}.nrep), 'LineWidth',2)
xlabel('Pitch (Hz)')
plot(F0s, rates_vel);
r = corrcoef(rates_vel, data_NT.rate);
R = r(1,2);
title([extractBefore(params_NT{1}.filename{1}, '.') ', R=' num2str(R)])
set(gca, 'fontsize', 16)
set(gca, 'xscale', 'log')
xticks([60 110 220 440])


%% Plot all as heatmap?



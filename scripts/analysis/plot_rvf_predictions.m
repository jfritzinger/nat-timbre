%%
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
hasRVF = cellfun(@(s) contains(s, 'R'), sessions.RVF);
hasBassoon = cellfun(@(s) contains(s, 'R'), sessions.Bassoon);
num_RVF = sum(hasRVF&hasBassoon);
RVFind = find(hasRVF&hasBassoon);

for ii = 1:num_RVF

	% Example session
	%putative = 'R29_TT3_P5_N02';
	putative = sessions.Putative_Units{RVFind(ii)};
	filename = sprintf('%s.mat', putative);
	load(fullfile(datapath,'neural_data', filename)), 'data';
	CF = sessions.CF(RVFind(ii));
	MTF_shape = sessions.MTF{RVFind(ii)};

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
	for ifile = 1:nfiles

		% Get single stimulus
		target = extractBefore(files{ifile}, '.');
		target_F0 = extractBetween(files{ifile}, 'ff.', '.wav');
		%savefile = sprintf('%s_F0_%s_Velocity.mat', target, target_F0{1});
		%load(fullfile(savepath, savefile),'vel_all', 'CFs_all')

		savefile = sprintf('%s_F0_%d_Velocity.mat', target, round(F0s(ifile)));
		load(fullfile(savepath, savefile), 'v', 'v_harms')
		vel_all = v;
		CFs_all = v_harms;

		% Find
		[~, CF_ind] = min(abs(CF-CFs_all));
		vel_CF(ifile) = vel_all(CF_ind);

	end

	%% Get rate / interpolated rate from RVF

	% Plot RVF
	params_RVF = data{5, 1};
	[fig, data_RVF] = plotPhysRVF([], params_RVF, []);
	close all

	%% Rearrange
	[sorted_vel, ind] = sort(data_RVF.velocities);
	sorted_rate = data_RVF.rate(ind);

	% Interpolate (linear)
	RVF_vel = -9:0.1:9;
	RVF_rates = interp1(sorted_vel, sorted_rate, RVF_vel, 'linear');

	% Match each velocity to a rate
	for ifiles = 1:nfiles

		% Find rate most similar to velocity
		this_vel = vel_CF(ifiles);
		[~, index] = min(abs(this_vel-RVF_vel));
		if index == 1
			rates_vel(ifiles) = RVF_rates(index+1);
		elseif index == length(RVF_vel)
			rates_vel(ifiles) = RVF_rates(index-1);
		else
			rates_vel(ifiles) = RVF_rates(index);
		end
	end

	% Plot both
	if ~isempty(data_NT.rate)
		r = corrcoef(rates_vel, data_NT.rate);
		R = r(1,2);
		R_all(ii) = R;
		R2_all(ii) = R^2;

		% Create matrices
		data_mat(ii,:) = data_NT.rate;
		model_mat(ii,:) = rates_vel;
	end
	

	% Plot natural timbre
	% figure
	% nexttile
	% hold on
	% errorbar(F0s,data_NT.rate, data_NT.rate_std/sqrt(params_NT{1}.nrep), 'LineWidth',2)
	% xlabel('Pitch (Hz)')
	% plot(F0s, rates_vel);
	% title([extractBefore(params_NT{1}.filename{1}, '.') ', R=' num2str(R)])
	% set(gca, 'fontsize', 16)

end


%% All together now 

figure
tiledlayout(1, 2)

nexttile
edges = linspace(-1, 1, 21);
histogram(R_all, edges)
xlim([-1 1])
xlabel('Correlation')
title('Correlation')
set(gca, 'fontsize', 16)

nexttile
edges = linspace(0, 1, 21);
histogram(R2_all, edges)
xlim([0 1])
xlabel('Variance Explained')
title('Variance Explained')
set(gca, 'fontsize', 16)

%% Find significant predictions 

shuffles = 1000;
alpha = 0.05;
[~, significant, correlation] = test_correlation_significance(model_mat, ...
	data_mat, shuffles, alpha);
sig_predictions = sum(significant);
fprintf('Sig predictions = %d/%d\n', sig_predictions, num_RVF);


%% FUNCTIONS 

function [p_values, is_significant, true_correlations] = ...
	test_correlation_significance(predictions,...
	actual_data, num_shuffles, alpha)
    % Input:
    % predictions: matrix of model predictions (150 neurons x num_samples)
    % actual_data: matrix of actual neural responses (150 neurons x num_samples)
    % num_shuffles: number of shuffles for permutation test (e.g., 1000)
    % alpha: significance level (e.g., 0.05)
	% From Claude AI, permutation testing
    
    num_neurons = size(predictions, 1);
    true_correlations = zeros(num_neurons, 1);
    null_correlations = zeros(num_neurons, num_shuffles);
    
    % Calculate true correlations
    for i = 1:num_neurons
        true_correlations(i) = corr(predictions(i,:)', actual_data(i,:)');
    end
    
    % Generate null distribution through shuffling
    for n = 1:num_neurons
        for s = 1:num_shuffles
            shuffled_data = actual_data(n,randperm(size(actual_data,2)));
            null_correlations(n,s) = corr(predictions(n,:)', shuffled_data');
        end
    end
    
    % Calculate p-values
    p_values = zeros(num_neurons, 1);
    for n = 1:num_neurons
        p_values(n) = mean(abs(null_correlations(n,:)) >= abs(true_correlations(n)));
    end
    
    % Determine significance
    is_significant = p_values < alpha;
end
%% calculate_null_model_distributions
clear  

%% Load in model and data and arrange into matrices
 
% Load in spreadsheet
addpath('/Users/jfritzinger/Projects/synth-timbre/scripts/helper-functions', '-end')
addpath '/Users/jfritzinger/Projects/nat-timbre/scripts/helper-functions'

[~, datapath, ~, ~] = getPaths();
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(datapath, 'data-cleaning', spreadsheet_name), 'PreserveVariableNames',true);
modelpath = '/Volumes/Nat-Timbre/data/manuscript';

% Find sessions for target synthetic timbre response
bin200(:,1) = cellfun(@(s) contains(s, 'R'), sessions.Bassoon);
bin200(:,2) = cellfun(@(s) contains(s, 'R'), sessions.Oboe);
isMTF = strcmp(sessions.MTF, 'BE')|strcmp(sessions.MTF, 'BS');
bin200_MTF = bin200 & isMTF;

has_data = bin200_MTF(:,1);
indices = find(has_data);
num_index = length(indices);
CFs = sessions.CF(indices);
for isesh = 1:num_index

	% Load in data
	putative = sessions.Putative_Units{indices(isesh)};
	load(fullfile(modelpath,'SFIE_model', [putative '_SFIE.mat']), 'SFIE', 'params_NT')
	load(fullfile(modelpath,'energy_model', [putative '_Energy.mat']), 'energy')
	load(fullfile(modelpath,'lat_inh_model', [putative '_Lat_Inh.mat']), 'lat_inh')
	load(fullfile(datapath, 'neural_data', [putative '.mat']))

	% Bassoon matrices
	if ~isempty(SFIE{2})

		params_NT = data{14, 2};
		data_NT = analyzeNT(params_NT);
		voice = data_NT.pitch_num<260 & data_NT.pitch_num>80;
		high = data_NT.pitch_num>260;

		SFIE_mat(isesh,:) = SFIE{2}.rate(high);
		energy_mat(isesh,:) = energy{2}.rate(high);
		lat_mat(isesh,:) = lat_inh{2}.rate(high);

		% Data
		data_mat(isesh, :) = data_NT.rate(high);
	end
	fprintf('%s done, %d percent done\n', putative, round(isesh/num_index*100))
end


%% Save matrices for ease of use later 

% Find and get rid of rows with 0s 
zero_rows = SFIE_mat(:,1)==0;
SFIE_mat(zero_rows,:) = [];
energy_mat(zero_rows,:) = [];
lat_mat(zero_rows,:) = [];
data_mat(zero_rows,:) = [];
CFs(zero_rows) = [];

% Save 
savename = 'Bassoon_High_Matrices.mat';
save(fullfile(modelpath, savename), 'SFIE_mat', 'energy_mat', 'lat_mat',...
	'data_mat', 'CFs')

% Load 
% savename = 'Bassoon_Voice_Matrices.mat';
% load(fullfile(modelpath, savename), 'SFIE_mat', 'energy_mat', 'lat_mat',...
% 	'data_mat', 'CFs')

%% Permutation testing to get significant predictions for SFIE 

shuffles = 1000;
alpha = 0.05;
[~, significant, correlation] = test_correlation_significance(SFIE_mat, ...
	data_mat, shuffles, alpha);

% Get significant predictions  
SFIE_sig = SFIE_mat(significant,:);
data_sig = data_mat(significant, :);
CFs_sig = CFs(significant);

% Plot results
figure('Position',[208,79,1321,774])
tiledlayout(4, 5)
nexttile
histogram(correlation)
title(['Significant: ' num2str(sum(significant)) '/' num2str(length(correlation))])

for ii = 1:sum(significant)
	nexttile
	plot(data_sig(ii,:));
	hold on
	plot(SFIE_sig(ii,:));
	correlation2 = corrcoef(data_sig(ii,:), SFIE_sig(ii,:));
	title(['R=' num2str(correlation2(1,2)) ', CF=' num2str(round(CFs_sig(ii)))])
end

% Plot scatter 
nexttile
hold on
scatter(CFs(significant), correlation(significant), 'filled',...
	'markeredgecolor', 'k')
scatter(CFs(~significant), correlation(~significant), 'filled',...
	'markeredgecolor', 'k')
set(gca, 'xscale', 'log')


%% Permutation testing to get significant predictions for energy 

shuffles = 1000;
alpha = 0.05;
[~, significant, correlation] = test_correlation_significance(energy_mat, ...
	data_mat, shuffles, alpha);

% Get significant predictions  
energy_sig = energy_mat(significant,:);
data_sig = data_mat(significant, :);
CFs_sig = CFs(significant);

% Plot results
figure('Position',[208,79,1321,774])
tiledlayout(4, 5)
nexttile
histogram(correlation)
title(['Significant: ' num2str(sum(significant)) '/' num2str(length(correlation))])

for ii = 1:sum(significant)
	nexttile
	plot(data_sig(ii,:));
	hold on
	plot(energy_sig(ii,:) .* (max(data_sig(ii,:))/max(energy_sig(ii,:))));
	correlation2 = corrcoef(data_sig(ii,:), energy_sig(ii,:));
	title(['R=' num2str(correlation2(1,2)) ', CF=' num2str(round(CFs_sig(ii)))])
end

% Plot scatter 
nexttile
hold on
scatter(CFs(significant), correlation(significant), 'filled',...
	'markeredgecolor', 'k')
scatter(CFs(~significant), correlation(~significant), 'filled',...
	'markeredgecolor', 'k')
set(gca, 'xscale', 'log')

%% Permutation testing to get significant predictions for lat inh 

shuffles = 1000;
alpha = 0.05;
[~, significant, correlation] = test_correlation_significance(lat_mat, ...
	data_mat, shuffles, alpha);

% Get significant predictions  
lat_sig = lat_mat(significant,:);
data_sig = data_mat(significant, :);
CFs_sig = CFs(significant);

% Plot results
figure('Position',[208,79,1321,774])
tiledlayout(4, 4)
nexttile
histogram(correlation)
title(['Significant: ' num2str(sum(significant)) '/' num2str(length(correlation))])

for ii = 1:sum(significant)
	nexttile
	plot(data_sig(ii,:));
	hold on
	plot(lat_sig(ii,:));
	correlation2 = corrcoef(data_sig(ii,:), lat_sig(ii,:));
	title(['R=' num2str(correlation2(1,2)) ', CF=' num2str(round(CFs_sig(ii)))])
end

% Plot scatter 
nexttile
hold on
scatter(CFs(significant), correlation(significant), 'filled',...
	'markeredgecolor', 'k')
scatter(CFs(~significant), correlation(~significant), 'filled',...
	'markeredgecolor', 'k')
set(gca, 'xscale', 'log')

%%

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
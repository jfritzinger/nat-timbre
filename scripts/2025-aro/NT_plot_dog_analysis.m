%% plot_dog_analysis

%% plot_dog_vs_gaussian
clear 

%% Load 

datapath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/';
load(fullfile(datapath, 'dog_analysis.mat'), "dog_analysis")


%% Do normal R^2 calculation

num_sesh = length(dog_analysis);
for ii = 1:num_sesh

	r = corrcoef(dog_analysis(ii).rate, dog_analysis(ii).dog_predicted);
	r2(ii) = r(1,2)^2;
	r_all(ii) = r(1,2);

	data_mat(ii, :) = dog_analysis(ii).rate;
	model_mat(ii, :) = dog_analysis(ii).dog_predicted;
end


%% Plot 
fontsize = 20;

figure('Position',[220,571,400,220])
nexttile
edges = linspace(0, 1, 21);
histogram(r2, edges )
xlim([0 1])
hold on
ylabel('# Neurons')
title('R^2 Fits')
xlabel('Variance Explained (R^2)')

[base, datapath] = getPaths();
load(fullfile(datapath, 'dog_analysis.mat'), "R2_gauss_all", "R2_dog_all", "dog_analysis")
histogram(R2_dog_all, edges, 'DisplayStyle','stairs', 'EdgeColor',[0.4 0.4 0.4],'LineWidth',1.5)
legend('Instrumental', 'Synthetic')
set(gca, 'fontSize', fontsize)

%% Get best example 
fontsize = 18;

datapath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/';
load(fullfile(datapath, 'dog_analysis.mat'), "R2_gauss_all", "R2_dog_all", "dog_analysis")
[max_r2, index] = max(r2);
pitch = categorical(dog_analysis(index).pitch);

figure('Position',[100 100 800 200])
tiledlayout(1, 3, 'Padding','compact')
nexttile([1, 2])
hold on
bar(pitch, dog_analysis(index).rate);
errorbar(pitch, dog_analysis(index).rate, dog_analysis(index).rate_std./...
	sqrt(30), 'LineStyle','none', 'LineWidth',1.5, 'Color','k');
plot(pitch, dog_analysis(index).dog_predicted, 'k', 'LineWidth',2);
xlabel('F0 (Hz)')
ylabel('Avg. Rate (sp/s)')
title('Example Neuron Fit')
legend('Data', '', 'DoG Fit', 'Location','northwest')
set(gca, 'fontsize', fontsize)

ax = gca;
categories = ax.XTickLabel;
new_labels = categories;
new_labels(2:2:end) = {''};  % Replace with empty strings
ax.XTickLabel = new_labels;


% Parameters
% g_exc, g_inh, s_exc, s_inh,  CF_exc, CF_inh
DOGparams = dog_analysis(index).dog_params;
Fs = 100000;
f = linspace(0, Fs/2, 100000);

% Set Parameters
s_exc = DOGparams(1);
s_inh = DOGparams(2);
sigma_exc = 10^DOGparams(3);
sigma_inh = 10^DOGparams(4);
CF_exc = 10^DOGparams(5);
CF_inh = 10^DOGparams(6);
gauss_exc = normpdf(f, CF_exc, sigma_exc);
gauss_inh = normpdf(f, CF_inh, sigma_inh);
gauss_exc = s_exc*(gauss_exc./max(gauss_exc));
gauss_inh = s_inh*(gauss_inh./max(gauss_inh));
W = gauss_exc - gauss_inh;

% % Plot to test
nexttile
%plot(gauss_exc)
hold on
%plot(-1*gauss_inh)
plot(W)
xlim([0 10000])
ticks = yticks;
yticklabels([0 1 2 3])
title('DoG Filter')
set(gca, 'fontsize', fontsize)
ylabel('Amplitude')
xlabel('Frequency (Hz)')
xline(dog_analysis(index).CF, '--', 'linewidth', 2)
set(gca, 'xscale', 'log')
xlim([300 10000])
legend('Filter', 'CF', 'Location','northwest')
xticks([200 500 1000 2000 5000 10000])
xticklabels([200 500 1000 2000 5000 10000]./1000)

%% Test significance 

shuffles = 1000;
alpha = 0.05;
[~, significant, correlation] = test_correlation_significance(model_mat, ...
	data_mat, shuffles, alpha);
num_sig = sum(significant);

% Get significant predictions  
lat_sig = model_mat(significant,:);
data_sig = data_mat(significant, :);
CFs_sig = CFs(significant);

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
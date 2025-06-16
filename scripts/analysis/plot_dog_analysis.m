%% plot_dog_analysis
clear 

%% Load 

[base, datapath] = getPathsNT();
load(fullfile(base, 'dog_analysis2.mat'), "dog_analysis")


%% Do normal R^2 calculation

num_sesh = length(dog_analysis);
for ii = 1:num_sesh
	r = corrcoef(dog_analysis(ii).rate, dog_analysis(ii).dog_predicted);
	r2(ii) = r(1,2)^2;
	r_all(ii) = r(1,2);
	data_mat(ii, :) = dog_analysis(ii).rate;
	model_mat(ii, :) = dog_analysis(ii).dog_predicted;

	r_g = corrcoef(dog_analysis(ii).rate, dog_analysis(ii).gaus_predicted);
	r2_g(ii) = r_g(1,2)^2;
end

%% Plot 
fontsize = 16;

figure('Position',[220,571,1153,230])
tiledlayout(1, 4, 'Padding','compact')
nexttile
edges = linspace(0, 1, 21);
histogram(r2, edges )
xlim([0 1])
hold on
ylabel('# Neurons')
title('R^2 Fits DoG')
xlabel('Variance Explained (R^2)')
set(gca, 'fontSize', fontsize)


nexttile
edges = linspace(0, 1, 21);
histogram(r2_g, edges)
xlim([0 1])
hold on
ylabel('# Neurons')
title('R^2 Fits Gaussian')
xlabel('Variance Explained (R^2)')
set(gca, 'fontSize', fontsize)

%% 
% Get best example 
[max_r2, index] = max(r2);
pitch = categorical(dog_analysis(index).pitch);
nexttile([1, 2])
hold on
bar(pitch, dog_analysis(index).rate);
errorbar(pitch, dog_analysis(index).rate, dog_analysis(index).rate_std./...
	sqrt(30), 'LineStyle','none', 'LineWidth',1.5, 'Color','k');
plot(pitch, dog_analysis(index).dog_predicted, 'k', 'LineWidth',2);
xlabel('F0 (Hz)')
ylabel('Avg. Rate (sp/s)')
title('Example Neuron Fit')
legend('Data', '', 'DoG Fit')
set(gca, 'fontsize', fontsize)

ax = gca;
categories = ax.XTickLabel;
new_labels = categories;
new_labels(2:2:end) = {''};  % Replace with empty strings
ax.XTickLabel = new_labels;

% Parameters 
nexttile
title('DoG Filter')
set(gca, 'fontsize', fontsize)
ylabel('Amplitude')
xlabel('Frequency (Hz)')


%% Test significance 

shuffles = 1000;
alpha = 0.05;
[~, significant, correlation] = test_correlation_significance(model_mat, ...
	data_mat, shuffles, alpha);

% Get significant predictions  
lat_sig = model_mat(significant,:);
data_sig = data_mat(significant, :);

%% Get all filters for fits above 0.4 
% 
% igood = find(r2>0.4);
% num_good = length(igood);
% for jj = 1:num_good
% 	index = igood(jj);
% 
% 	if ~isempty(dog_analysis(index).dog_params)
% 		figure('Position',[560,615,956,233])
% 		tiledlayout(1, 3)
% 		pitch = categorical(dog_analysis(index).pitch);
% 		nexttile([1, 2])
% 		hold on
% 		bar(pitch, dog_analysis(index).rate);
% 		errorbar(pitch, dog_analysis(index).rate, dog_analysis(index).rate_std./...
% 			sqrt(30), 'LineStyle','none', 'LineWidth',1.5, 'Color','k');
% 		plot(pitch, dog_analysis(index).dog_predicted, 'k', 'LineWidth',2);
% 		xlabel('F0 (Hz)')
% 		ylabel('Avg. Rate (sp/s)')
% 		title(['Example Neuron Fit, R^2=' num2str(round(r2(index), 2))])
% 		legend('Data', '', 'DoG Fit')
% 		set(gca, 'fontsize', fontsize)
% 
% 		ax = gca;
% 		categories = ax.XTickLabel;
% 		new_labels = categories;
% 		new_labels(2:2:end) = {''};  % Replace with empty strings
% 		ax.XTickLabel = new_labels;
% 
% 		% Get variables and plot
% 		% g_exc, g_inh, s_exc, s_inh,  CF_exc, CF_inh
% 		DOGparams = dog_analysis(index).dog_params;
% 		Fs = 100000;
% 		f = linspace(0, Fs/2, 100000);
% 
% 		% Set Parameters
% 		s_exc = DOGparams(1);
% 		s_inh = DOGparams(2);
% 		sigma_exc = 10^DOGparams(3);
% 		sigma_inh = 10^DOGparams(4);
% 		CF_exc = 10^DOGparams(5);
% 		%CF_inh = 10^DOGparams(6);
% 		gauss_exc = normpdf(f, CF_exc, sigma_exc);
% 		gauss_inh = normpdf(f, CF_exc, sigma_inh);
% 		gauss_exc = s_exc*(gauss_exc./max(gauss_exc));
% 		gauss_inh = s_inh*(gauss_inh./max(gauss_inh));
% 		W = gauss_exc - gauss_inh;
% 
% 		% Plot to test
% 		nexttile
% 		plot(gauss_exc)
% 		hold on
% 		plot(-1*gauss_inh)
% 		plot(W)
% 		xlim([0 10000])
% 		title('DoG Filter')
% 		set(gca, 'fontsize', fontsize)
% 		ylabel('Amplitude')
% 		xlabel('Frequency (Hz)')
% 	end
% end

%% Plot Gaussian vs DoG 

sig_ind = [dog_analysis.p_value] < 0.05;
non_ind = [dog_analysis.p_value] >= 0.05;

figure
R2_dog = abs([dog_analysis.R2_dog]);
R2_gaus = abs([dog_analysis.R2_gauss]);

scatter(R2_gaus(sig_ind), R2_dog(sig_ind), 'filled', 'MarkerFaceAlpha',0.5, 'MarkerEdgeColor','k')
hold on
scatter(R2_gaus(non_ind), R2_dog(non_ind), 'filled', 'MarkerFaceAlpha',0.5, 'MarkerEdgeColor','k')
xlim([0 1])
ylim([0 1])
plot([0 1], [0 1], 'k')
xlabel('R^2 Gaussian')
ylabel('R^2 DoG')
legend('Sig.', 'Non sig.')

%% Plot scatter plots of fit parameters, similar to WB-TIN and synth timbre 
figure 
scattersize = 30;
legsize = 14;


% Get good fits
good_fit = [dog_analysis.R2_dog]>0.4;

% Plot 
all_dog_params = [dog_analysis(good_fit).dog_params];
all_dog_params = reshape(all_dog_params, 5,[])'; % 6 for OG, 5 for new
CFs = [dog_analysis(good_fit).CF];

% Un-log CF_exc, CF_inh
CF_exc = 10.^all_dog_params(:,5);
s_exc = 10.^all_dog_params(:,3);
s_inh = 10.^all_dog_params(:,4);
g_exc = all_dog_params(:,1);
g_inh = all_dog_params(:,2);

% Scatter plot of ratio of inhibition to excitation sigma and strengths 
ratio_sigma = log10(s_inh./s_exc);
ratio_g = log10(g_inh./g_exc);
nexttile
hold on
scatter(ratio_sigma, ratio_g, scattersize, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.5)
xline(0)
yline(0)
xlabel('Log BW Ratio (\sigma_i_n_h/\sigma_e_x_c)')
ylabel('Log Str Ratio (g_i_n_h/g_e_x_c)')
grid on
xlim([-2.1 2.3])
ylim([-1.5 1.5])
set(gca, 'fontsize', fontsize)

% Number of units in each quadrant
q1 = sum(ratio_sigma>0 & ratio_g>0);
q2 = sum(ratio_sigma<0 & ratio_g>0);
q3 = sum(ratio_sigma<0 & ratio_g<0);
q4 = sum(ratio_sigma>0 & ratio_g<0);

xL=xlim(gca);
yL=ylim(gca);
text(gca, 0.95*xL(1),0.99*yL(2),sprintf('n=%d', q2),...
	'HorizontalAlignment','left','VerticalAlignment','top', 'FontSize',legsize)
text(gca, 0.95*xL(2),0.99*yL(2),sprintf('n=%d', q1),...
	'HorizontalAlignment','right','VerticalAlignment','top', 'FontSize',legsize)
text(gca, 0.95*xL(1),0.99*yL(1),sprintf('n=%d', q3),...
	'HorizontalAlignment','left','VerticalAlignment','bottom', 'FontSize',legsize)
text(gca, 0.95*xL(2),0.99*yL(1),sprintf('n=%d', q4),...
	'HorizontalAlignment','right','VerticalAlignment','bottom', 'FontSize',legsize)



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




%% plot_dog_analysis

%% plot_dog_vs_gaussian
clear 

%% Load 

datapath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/';
load(fullfile(datapath, 'dog_analysis.mat'), "R2_gauss_all", "R2_dog_all", "dog_analysis")


%% Do normal R^2 calculation

num_sesh = length(dog_analysis);
for ii = 1:num_sesh

	r = corrcoef(dog_analysis(ii).rate, dog_analysis(ii).dog_predicted);
	r2(ii) = r(1,2)^2;
	r_all(ii) = r(1,2);

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
title('R^2 Fits')
xlabel('Variance Explained (R^2)')

[base, datapath] = getPaths();
load(fullfile(datapath, 'dog_analysis.mat'), "R2_gauss_all", "R2_dog_all", "dog_analysis")
histogram(R2_dog_all, edges, 'DisplayStyle','stairs', 'EdgeColor',[0.4 0.4 0.4],'LineWidth',1.5)
legend('Instrumental Timbre', 'Synthetic Timbre')
set(gca, 'fontSize', fontsize)

% Get best example 
datapath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/';
load(fullfile(datapath, 'dog_analysis.mat'), "R2_gauss_all", "R2_dog_all", "dog_analysis")
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


%%

figure
edges = linspace(-1, 1, 21);
histogram(r_all, edges )
xlim([-1 1])
hold on
ylabel('# Neurons')
title('Correlation')
set(gca, 'fontSize', 16)


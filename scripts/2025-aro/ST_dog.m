%% ST_dog%% dog_analysis_plots 
clear 

%% Load in fits  

[base, datapath] = getPaths();
load(fullfile(datapath, 'dog_analysis.mat'), "R2_gauss_all", "R2_dog_all", "dog_analysis")


%% Plot example fit 

figure('Position',[53,540,450,500])
tiledlayout(2, 1)
fontsize = 20;

% Example: R24_TT2_P13_N02, CF = 1150Hz, BS
putative = 'R24_TT2_P13_N02';
ind = cellfun(@(d) strcmp(d, putative), {dog_analysis.putative});

% Load in to get spont rate
load(fullfile(datapath, 'neural_data', [putative '.mat']))
params_RM = data{2, 2};
data_RM = analyzeRM(params_RM);
spont = data_RM.spont;

% Plot 
nexttile
hold on
plot(dog_analysis(ind).fpeaks, dog_analysis(ind).rate, 'k','LineWidth',2);
plot(dog_analysis(ind).fpeaks, dog_analysis(ind).dog_predicted, ...
	'LineWidth',2, 'color', 'b');
plot(dog_analysis(ind).fpeaks, dog_analysis(ind).gaus_predicted,...
	'LineWidth',2, 'color', '#1b9e77');
xline(dog_analysis(ind).CF, '--', 'LineWidth',2)
yline(spont, 'k')
ylabel('Avg. Rate (sp/s)')
xlabel('Spectral Peak Freq. (Hz)')
title('Example Fit')
set(gca, 'fontSize', fontsize)
hleg = legend('Data', 'DoG', 'Gaussian', 'CF', 'fontsize', 16);
hleg.ItemTokenSize = [16, 6];
grid on
xlim([dog_analysis(ind).fpeaks(1) dog_analysis(ind).fpeaks(end)])

%% Plot adjusted R^2 values 

sig = [dog_analysis.p_value]<0.05;
notsig = [dog_analysis.p_value]>0.05;

nexttile
scatter(R2_gauss_all(sig), R2_dog_all(sig),60, 'filled',...
	'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.6, 'MarkerFaceColor','b')
hold on
scatter(R2_gauss_all(notsig), R2_dog_all(notsig), 60, 'filled',...
	'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.6, 'MarkerFaceColor',[0.4 0.4 0.4])
xlim([0 1])
ylim([0 1])
xticks([0 0.2 0.4 0.6 0.8 1])
yticks(0:0.2:1)
grid on
plot([0 1], [0 1], 'k')
xlabel('Gaussian Adjusted R^2')
ylabel('DoG Adjusted R^2')
title('DoG vs Gaussian Comparison')
set(gca, 'fontSize', fontsize)
msg = sprintf('%d sig.', sum(sig));
msg2 = sprintf('%d not sig.', sum(notsig));
legend(msg, msg2, 'Location','southeast')

%% Plot F-test for MSE of fits 

% num_sesh = length(dog_analysis);
% p_all = NaN(num_sesh, 1);
% CF_all = NaN(num_sesh, 1);
% for ii = 1:num_sesh
% 
% 	current_dog = dog_analysis(ii).dog_predicted;
% 	current_gauss = dog_analysis(ii).gaus_predicted;
% 	rate = dog_analysis(ii).rate;
% 	CF_all(ii) = dog_analysis(ii).CF;
% 
% 	% Comparing based on how close the curves are to data 
% 	p_value = ftest(rate, current_gauss, current_dog);
% 	p_all(ii) = log(p_value);
% 
% end
% 
% sig = p_all<log(0.05);
% not_sig = p_all>log(0.05);
% 
% nexttile
% scatter(CF_all(sig), p_all(sig), 'filled', "MarkerEdgeColor",'k')
% hold on
% scatter(CF_all(not_sig), p_all(not_sig), 'filled', "MarkerEdgeColor",'k')
% yline(log(0.05), '--', 'linewidth', 2.5)
% set(gca, 'XScale', 'log')
% ylabel('Log p-value from F-test')
% xlabel('CF')
% set(gca, 'fontsize', fontsize)
% title({'Is DoG significantly', 'better than Gaussian?'})
% msg = sprintf('%d sig.', sum(sig));
% msg2 = sprintf('%d not sig.', sum(not_sig));
% legend(msg, msg2, 'p=0.05')

% Annotations 
% text(0.65, 0.15, msg, 'Units', 'normalized', ...
% 	'VerticalAlignment', 'top', 'FontSize',16)
% text(0.65, 0.09, msg2, 'Units', 'normalized', ...
% 	'VerticalAlignment', 'top', 'FontSize',16)

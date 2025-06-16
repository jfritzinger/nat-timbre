%% plot_neuron_time
clear

%% Load in data
target = 'Bassoon';

pitch = getF0s(target);
msg{1} = sprintf('%0.0f - %0.0f Hz F0', pitch(1), pitch(20));
msg{2} = sprintf('%0.0f - %0.0f Hz F0', pitch(21), pitch(end));

[base, ~, ~, ~] = getPathsNT();
load(fullfile(base, 'model_comparisons', ['Neuron_Time_F0_' target '.mat']), ...
	"neuron_time_F0")
% load(fullfile(base, 'model_comparisons', 'Neuron_Time_Timbre_All.mat'),...
% 	"neuron_time_timbre")
% neuron_time_F0 = neuron_time_timbre;

%% Plot accuracy of each neuron
figure('Position',[560,618,798,230])
tiledlayout(1, 3)
linewidth = 1;

accuracy(1,:) = [neuron_time_F0.accuracy]*100;
accuracy(2,:) = [neuron_time_F0.accuracy_low]*100;
accuracy(3,:) = [neuron_time_F0.accuracy_high]*100;

% Plot overall histogram
nexttile
edges = linspace(0, 100, 30);
histogram(accuracy(1,:),edges)
hold on
chance = 1/size(neuron_time_F0(1).C, 1)*100;
xline(chance, 'k', 'LineWidth',linewidth)
xline(mean(accuracy(1,:)), 'r', 'LineWidth',2)
scatter(max(accuracy(1,:)), 0, 'filled')
ylabel('# Neurons')
xlabel('Prediction Accuracy (%)')
title([target ' Prediction of F0'])
xlim([0 100])
legend('', sprintf('Chance = %.2f%%', chance), ...
	sprintf('Mean = %.2f%%', mean(accuracy(1,:))), ...
	sprintf('Best = %.2f%%', max(accuracy(1,:))))

% Split histogram into two sections
nexttile
edges = linspace(0, 100, 30);
histogram(accuracy(2,:),edges, 'FaceAlpha',0.5, 'FaceColor','b')
hold on
histogram(accuracy(3,:),edges, 'FaceAlpha',0.5, 'FaceColor','r')
xline(chance, 'k', 'LineWidth',linewidth)
xline(mean(accuracy(2,:)), 'b', 'LineWidth',1)
xline(mean(accuracy(3,:)), 'r', 'LineWidth',1)
ylabel('# Neurons')
xlabel('Prediction Accuracy (%)')
title('Low vs high F0 Predictions')
xlim([0 100])
legend(msg{1}, msg{2}, ...
	sprintf('Chance = %.2f%%', chance), ...
	sprintf('Mean = %.2f%%', mean(accuracy(2,:))), ...
	sprintf('Mean = %.2f%%', mean(accuracy(3,:))))


% Find indices of the best neurons 
[temp,originalpos] = sort(accuracy(1,:), 'descend' );
n = temp(1:3);
best_ind=originalpos(1:3);

% Plot confusion matrix of best neuron
putatives = neuron_time_F0(best_ind(1)).putative;
nexttile
confusionchart(neuron_time_F0(best_ind(1)).C)
title(sprintf('Best prediction, %s, %0.02f%%', putatives, temp(1)))


%% Test VS vs accuracy 

load(fullfile(base,'model_comparisons', 'Data_NT_3.mat'), 'nat_data')

% Get subset of units 
VS_all = zeros(3,1);
for ii = 1:297
	VS_1 = nat_data(ii).bass_VS;
	if ~isempty(VS_1)
		VS_all(ii,:) = mean(VS_1);
	end
end
VS_all(VS_all==0) = [];

figure
tiledlayout(3, 3)

for ii = 1:3
	nexttile
	histogram(accuracy(ii,:))
	accuracy_all = accuracy(ii,:);

	nexttile
	histogram(VS_all)

	nexttile
	scatter(VS_all, accuracy(ii,:));
	hold on

	mdl = fitlm(VS_all, accuracy(ii,:));
	x = linspace(0, 0.6, 20);
	y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
	plot(x, y, 'r')
	hleg = legend('Neuron', ...
		sprintf('y = %0.2f*x+%0.2f, p=%0.04f', mdl.Coefficients{2, 1}, ...
		mdl.Coefficients{1, 1},mdl.Coefficients{2,4}));
	ylabel('Accuracy (%)')
	xlabel('Vector Strength')
end

%% Test accuracy vs neuron characteristics, like CF and MTF 

% Plot CF vs accuracy 
figure

nexttile
CFs = [neuron_time_F0.CF];
scatter(CFs, accuracy(1,:));
set(gca, 'xscale', 'log')
hold on
mdl = fitlm(CFs, accuracy(1,:));
x = linspace(300, 10000, 20);
y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
plot(x, y, 'r')
hleg = legend('Neuron', ...
	sprintf('y = %0.2f*x+%0.2f, p=%0.04f', mdl.Coefficients{2, 1}, ...
	mdl.Coefficients{1, 1},mdl.Coefficients{2,4}));

nexttile
scatter(CFs, accuracy(2,:));
set(gca, 'xscale', 'log')

nexttile
scatter(CFs, accuracy(2,:));
set(gca, 'xscale', 'log')


%% Plot MTF group vs accuracy 

MTFs = {neuron_time_F0.MTF};
MTF_types = unique(MTFs);

figure
hold on
for iMTF = 1:5
	ind = strcmp(MTFs, MTF_types{iMTF});
	accur = accuracy(1,ind);
	num_units = length(accur);

	swarmchart(ones(num_units, 1)*iMTF, accur)

	mean_vals(iMTF) = mean(accur);
	std_vals(iMTF) = std(accur)/sqrt(length(accur));
end
errorbar(1:5, mean_vals, std_vals, 'k')
xticks(1:5)
xticklabels(MTF_types)
xlabel('MTF Groups')
ylabel('Accuracy')

tableMTF = table(MTFs', accuracy(1,:)');
anova(tableMTF, 'Var2')
[~,~,stats] = anova1(accuracy(1,:), MTFs);
[c,~,~,gnames] = multcompare(stats);
grid on


%% Plot PC2 vs accuracy


VS_all = zeros(3,1);
for ii = 1:297
	VS_1 = nat_data(ii).bass_VS;
	if ~isempty(VS_1)
		PC2_score(ii,:) = nat_data(ii).RVF_PC2;
	else
		PC2_score(ii,:) = 0;
	end
end
PC2_score(PC2_score==0) = [];


% Plot PCA2 score vs beta weights 
figure
scatter(PC2_score, accuracy(1,:), 'filled', 'MarkerEdgeColor','k')
xlabel('PC2 RVF score')
ylabel('Accuracy')
hold on

mdl = fitlm(PC2_score, accuracy(1,:));
x = linspace(-2, 1, 20);
y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
plot(x, y, 'r')
hleg = legend('Neuron', ...
	sprintf('p=%0.04f',mdl.Coefficients{2,4}), 'fontsize', 7, ...
	'location', 'northwest', 'box', 'off');
hleg.ItemTokenSize = [8, 8];
grid on



%% Rate vs timing comparison 

% Get best single-unit timbre neurons 
[base, ~, ~, ppi] = getPathsNT();
filepath = fullfile(base, 'model_comparisons','Neuron_Rate_F0_Bassoon.mat');
load(filepath, "neuron_rate_F0")
accuracy_rate = [neuron_rate_F0.accuracy]*100;
accuracy_time = accuracy(1,:);

figure
hold on
scatter(accuracy_rate, accuracy_time, 'filled', 'MarkerEdgeColor','k')
xlabel('Rate Accuracy')
ylabel('Timing Accuracy')
plot([0 100], [0 100], 'k')

xlim([0 20])
ylim([0 60])

mdl = fitlm(accuracy_rate, accuracy_time);
x = linspace(0, 60, 10);
y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
plot(x, y, 'r')
hleg = legend('Neuron', 'p=7.1671e-25');
hleg.ItemTokenSize = [8, 8];
title('Rate vs Timing')

figure
edges = linspace(0, 20, 21);
histogram(accuracy_rate,edges)

figure
edges = linspace(0, 60, 21);
histogram(accuracy_time, edges)

%% Plot imagesc for diag for each neuron 

[~,best_ind] = sort(abs(accuracy(1,:)), 'ascend' );
for ii = 1:287
	C_acc(ii,:) = diag(neuron_time_F0(best_ind(ii)).C);
end

x = getF0s(target);
y = 1:287;

figure
pcolor(x, y, C_acc, 'EdgeColor','none')
set(gca, 'xscale', 'log')
xticks([60 100 200 350 550])
xlabel('F0 (Hz)')
ylabel('All Neurons')
colorbar
title('Correct predictions of\newlineeach F0 for all neurons')

%% Load and plot results for timbre all F0s together 
% 
% % Load in data 
% [base, datapath, savepath, ppi] = getPathsNT();
% filename = 'Neuron_Time_Timbre_All.mat';
% filepath = fullfile(base, 'model_comparisons',filename);
% load(filepath, "neuron_time_timbre")
% 
% % Figure
% figure('Position',[560,585,1023,263])
% tiledlayout(1, 2)
% linewidth = 2;
% 
% nexttile
% accuracy = [neuron_time_timbre.accuracy]*100;
% CFs = [neuron_time_timbre.CF];
% edges = linspace(0, 100, 101);
% histogram(accuracy,edges)
% hold on
% 
% % Plot chance line
% xline(50, 'k', 'LineWidth',linewidth)
% 
% % Mean and median
% xline(mean(accuracy), 'r', 'LineWidth',linewidth)
% xline(median(accuracy), 'r--', 'LineWidth',linewidth)
% 
% % Best
% scatter(max(accuracy), 0, 'filled')
% 
% % Labels
% legend('', sprintf('Chance = 50%%'), ...
% 	sprintf('Mean = %.2f%%', mean(accuracy)), ...
% 	sprintf('Median = %.2f%%', median(accuracy)),...
% 	sprintf('Best = %.2f%%', max(accuracy)))
% ylabel('# Neurons')
% xlabel('Prediction Accuracy (%)')
% title('Instrument Predictions')
% xlim([0 100])
% grid on
% 
% % Plot accuracy vs CF
% nexttile
% scatter(CFs, accuracy, 'filled')
% set(gca, 'xscale', 'log')
% xlabel('CFs (Hz)')
% ylabel('Accuracy')


%% FUNCTIONS 

function pitch = getF0s(target)

[base, ~, ~, ~] = getPathsNT();
tuning = readtable(fullfile(base, 'Tuning.xlsx')); % Load in tuning
listing = dir(fullfile(base, 'waveforms', '*.wav'));
target_WAV = arrayfun(@(n) contains(listing(n).name, target), 1:numel(listing), 'UniformOutput', false);
wav_nums =  find(cell2mat(target_WAV));
d = dir(fullfile(base,'waveforms', '*.wav'));
all_files = sort({d.name});
nfiles = length(wav_nums);

for i = 1:nfiles
	files{1,i} = all_files{wav_nums(i)};
end

% Sort by frequency of pitch
index = [];
note_names = extractBetween(files, 'ff.','.');
for ii = 1:nfiles % Find index of each note in tuning spreadsheet
	index(ii) = find(strcmp(note_names(ii), tuning.Note));
end
pitch_order = tuning.Frequency(index); % Get freqs of each note
[~, order] = sort(pitch_order); % Sort freqs
pitch = pitch_order(order);

end
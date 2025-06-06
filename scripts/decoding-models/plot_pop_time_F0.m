%% plot_pop_time_F0
clear

%% Set up figure

figure
tiledlayout(2, 3)

%% Load in data 
target = 'Bassoon';

base = getPathsNT();
load(fullfile(base, 'model_comparisons', 'pop_timing_F0_bassoon_subset.mat'), ...
	"accur_all","C_all", "num_neurons")

%target = 'Oboe';

if ismac
	fpath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data';
else
	fpath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\';
end
tuning = readtable(fullfile(fpath, 'Tuning.xlsx')); % Load in tuning

% Get bassoon stimulus
listing = dir(fullfile(fpath, 'waveforms', ['*' target '*.wav']));
files = {listing.name};
note_names = extractBetween(files, 'ff.', '.');
[~, index] = ismember(note_names, tuning.Note);
F0s1 = tuning.Frequency(index);
[F0s, order] = sort(F0s1);

%% Plot 

nmodels = length(num_neurons);
mean_acc = max(accur_all,[],2);
std_acc = std(accur_all, [], 2);


nexttile
errorbar(num_neurons(1:nmodels/2), mean_acc(1:nmodels/2), std_acc(1:nmodels/2)/sqrt(10));
hold on
errorbar(num_neurons(nmodels/2+1:end), mean_acc(nmodels/2+1:end), std_acc(nmodels/2+1:end)/sqrt(10))
xlabel('Number of Neurons in Model')
ylabel('Accuracy')
legend('Best', 'Worst')
title('Model Accuracy')
ylim([0 1])

%% 

C_diag = NaN(nmodels, 40);
for ii = 1:nmodels 
	C_diag(ii,:) = diag(C_all{ii,20});
end
y = num_neurons(1:12);
x = F0s;

nexttile
pcolor(x, y, C_diag(1:12,:), 'EdgeColor','none')
set(gca, 'xscale', 'log')
clim([0 20])
ylabel('# Neurons in Model')
xlabel('F0s')
title('Best Units')
a=colorbar;
a.Label.String = '# Accurate Predictions'; 
xticks([55 110 220 440])

nexttile
pcolor(x, y, C_diag(13:24,:), 'EdgeColor','none')
set(gca, 'xscale', 'log')
clim([0 20])
ylabel('# Neurons in Model')
xlabel('F0s')
title('Worst Units')
a=colorbar;
a.Label.String = '# Accurate Predictions'; 
xticks([55 110 220 440])

%% Annotate 




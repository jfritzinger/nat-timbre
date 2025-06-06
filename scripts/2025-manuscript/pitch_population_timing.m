%% plot_pop_time_F0
clear
save_fig = 0;


%% Set up figure

[base, ~, ~, ppi] = getPathsNT();
figure('Position',[50, 50, 6*ppi, 3*ppi])
linewidth = 1;
fontsize = 8;
legsize = 7;
labelsize = 12;
scattersize = 10;

%% Load in data

hind = [1 4];
for iinstru = 1:2

	target = 'Bassoon';
	%target = 'Oboe';
	load(fullfile(base, 'model_comparisons', 'pop_timing_F0_bassoon_subset.mat'), ...
		"accur_all","C_all", "num_neurons")

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
	mean_acc = mean(accur_all,2);
	std_acc = std(accur_all, [], 2);

	h(hind(iinstru)) = subplot(2, 3, hind(iinstru));
	errorbar(num_neurons(1:nmodels/2), mean_acc(1:nmodels/2), std_acc(1:nmodels/2)/sqrt(10));
	hold on
	errorbar(num_neurons(nmodels/2+1:end), mean_acc(nmodels/2+1:end), std_acc(nmodels/2+1:end)/sqrt(10))
	xlabel('# Neurons in Model')
	ylabel('Accuracy')
	hleg = legend('Best', 'Worst', 'fontsize', legsize);
	hleg.ItemTokenSize = [8, 8];
	title('Model Accuracy')
	ylim([0 1.2])
	grid on
	set(gca, 'FontSize', fontsize)
	box off

	%%

	C_diag_all = NaN(nmodels, 40);
	for ii = 1:nmodels

		for irep = 1:10
			C_diag(irep,:) = diag(C_all{ii,irep});
		end

		C_diag_all(ii,:) = mean(C_diag);
	end
	y = num_neurons(1:12);
	x = F0s;

	h(hind(iinstru)+1) = subplot(2, 3, hind(iinstru)+1);
	pcolor(x, y, C_diag_all(1:12,:), 'EdgeColor','none')
	set(gca, 'xscale', 'log')
	clim([0 20])
	ylabel('# Neurons in Model')
	xlabel('F0s')
	title('Best Units')
	a=colorbar;
	a.Label.String = '# Accurate Predictions';
	xticks([55 110 220 440])
	set(gca, 'FontSize', fontsize)
	box off

	%%

	h(hind(iinstru)+2) = subplot(2, 3, hind(iinstru)+2);
	pcolor(x, y, C_diag_all(13:24,:), 'EdgeColor','none')
	set(gca, 'xscale', 'log')
	clim([0 20])
	ylabel('# Neurons in Model')
	xlabel('F0s')
	title('Worst Units')
	a=colorbar;
	a.Label.String = '# Accurate Predictions';
	xticks([55 110 220 440])
	set(gca, 'FontSize', fontsize)
	box off
end

%% Arrange positions

left = [0.12 0.43 0.74]; 
bottom = linspace(0.13, 0.6, 2);
height = 0.3;
width = 0.15;

set(h(1), 'position', [left(1) bottom(2) 0.21 height])
set(h(2), 'position', [left(2) bottom(2) width height])
set(h(3), 'position', [left(3) bottom(2) width height])
set(h(4), 'position', [left(1) bottom(1) 0.21 height])
set(h(5), 'position', [left(2) bottom(1) width height])
set(h(6), 'position', [left(3) bottom(1) width height])


%% Annotate

annotation("textbox", [0.05517 0.6 0.1386 0.1088], "String", "Bassoon",...
	"FontSize", 12, "FontWeight", "bold", "EdgeColor", "none", "Rotation",90)
annotation("textbox", [0.0348 0.18 0.09778 0.05666], "String", "Oboe",...
	"FontSize", 12, "FontWeight", "bold", "EdgeColor", "none", "Rotation",90)

%% Save figure

if save_fig == 1
	filename = 'figsomething';
	save_figure(filename)
end
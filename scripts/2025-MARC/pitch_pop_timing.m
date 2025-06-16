%% plot_pop_time_F0
clear
save_fig = 1;


%% Set up figure

[base, ~, ~, ppi] = getPathsNT();
figure('Position',[50, 50, 10.3*ppi, 6.4*ppi])
linewidth = 2;
fontsize = 18;
titlesize = 20;
legsize = 20;
scattersize = 40;

%% Load in data

hind = [1 3];
targets = {'Bassoon', 'Oboe'}; 
for iinstru = 1:2

	target = targets{iinstru};
	if iinstru == 2
		load(fullfile(base, 'model_comparisons', 'pop_timing_F0_Oboe_subset2.mat'), ...
			"accur_all","C_all", "num_neurons", "mean_acc", "std_acc")
		accur_all2 = accur_all;
		C_all2 = C_all;
		num_neurons2 = num_neurons;
		mean_acc2 = mean_acc;
		std_acc2 = std_acc;
	end
	load(fullfile(base, 'model_comparisons', ['pop_timing_F0_' target '_subset.mat']), ...
		"accur_all","C_all", "num_neurons")
	load(fullfile(base, 'model_comparisons', ['Pop_Rate_F0_' target '3.mat']),...
		"pop_rate_F0")
	

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

	if iinstru == 1
	nmodels = length(num_neurons);
	mean_acc = mean(accur_all,2);
	std_acc = std(accur_all, [], 2);
	else
		nmodels = length(num_neurons) + length(num_neurons2);
		mean_acc1 = mean(accur_all,2);
		mean_acc = [mean_acc1(1:12); mean_acc2(1:3); mean_acc1(13:end); mean_acc2(4:6)];
		std_acc1 = std(accur_all, [], 2);
		std_acc = [std_acc1(1:12); std_acc2(1:3); std_acc1(13:end); std_acc2(4:6)];
		num_neurons = [num_neurons(1:12) num_neurons2(1:3) num_neurons(13:24) num_neurons2(4:6)];
	end

	h(hind(iinstru)) = subplot(2, 2, hind(iinstru));
	errorbar(num_neurons(1:nmodels/2), mean_acc(1:nmodels/2), std_acc(1:nmodels/2)/sqrt(10), 'LineWidth',linewidth);
	hold on
	%errorbar(num_neurons(nmodels/2+1:end), mean_acc(nmodels/2+1:end), std_acc(nmodels/2+1:end)/sqrt(10))
	xlabel('# Neurons in Model')
	ylabel('Accuracy')
	% hleg = legend('Best', 'Worst', 'fontsize', legsize);
	% hleg.ItemTokenSize = [8, 8];
	title([target ' Model Accuracy'])
	ylim([0 1])
	grid on
	set(gca, 'FontSize', fontsize)
	box off
	yticks(0:0.2:1)
	xticks(0:5:100)
	if iinstru == 2
		xticklabels({'1', '', '10', '', '20', '', '30', '', '40', '', '50', '', '60', '', '70', '', '80', '', '90', '', '100'})
	end

	%%

	C_diag_all = NaN(nmodels, length(F0s));
	if iinstru == 1
		for ii = 1:nmodels
			C_diag = NaN(10, length(F0s));
			for irep = 1:10
				C_diag(irep,:) = diag(C_all{ii,irep});
			end

			C_diag_all(ii,:) = mean(C_diag);
		end
	else
		iii = 1;
		iiii = 1;
		for ii = 1:nmodels
			if ~ismember(ii, [13, 14, 15, 28, 29, 30])
				C_diag = NaN(10, length(F0s));
				for irep = 1:10
					C_diag(irep,:) = diag(C_all{iiii,irep});
				end
				iiii = iiii +1;
				C_diag_all(ii,:) = mean(C_diag);
			else
				C_diag_all(ii,:) = diag(C_all2{iii});
				iii = iii +1;
			end
		end
	end
	y = F0s;
	x = num_neurons(1:nmodels/2);

	h(hind(iinstru)+1) = subplot(2, 2, hind(iinstru)+1);
	pcolor(x, y, C_diag_all(1:nmodels/2,:)', 'EdgeColor','none')
	set(gca, 'yscale', 'log')
	clim([0 20])
	xlabel('# Neurons in Model')
	%shading interp
	ylabel('F0 (kHz)')
	title(['Accuracy vs F0'])
	a=colorbar;
	a.Label.String = '# Accurate Predictions';
	yticks([50 100 250 500 1000 1500])
	yticklabels([50 100 250 500 1000 1500]/1000)
	set(gca, 'FontSize', fontsize)
	box off
	%yticks(0:0.2:1)
	%xticks(0:5:40)

end

%% Arrange positions

left = [0.12 0.55]; 
bottom = linspace(0.13, 0.62, 2);
height = 0.3;
width = 0.3;

set(h(1), 'position', [left(1) bottom(2) width height])
set(h(2), 'position', [left(2) bottom(2) width height])
set(h(3), 'position', [left(1) bottom(1) width height])
set(h(4), 'position', [left(2) bottom(1) width height])

%% Save figure

if save_fig == 1
	filename = 'pitch_pop_timing';
	save_figure_MARC(filename)
end
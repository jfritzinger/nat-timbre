%% plot_pop_time_F0
clear
save_fig = 0;


%% Set up figure

[base, ~, ~, ppi] = getPathsNT();
figure('Position',[50, 50, 4*ppi, 3.7*ppi])
linewidth = 1;
fontsize = 8;
legsize = 7;
labelsize = 12;
scattersize = 10;

%% Load in data
tuning = readtable(fullfile(base, 'Tuning.xlsx')); % Load in tuning
hind = [1 3 5];
targets = {'Bassoon', 'Oboe', 'Invariant'}; 
for iinstru = 1:3

	target = targets{iinstru};
	if iinstru == 1
		% load(fullfile(base, 'model_comparisons', ['pop_timing_F0_' target '_subset.mat']), ...
		% 	"accur_all","C_all", "num_neurons")
		% load(fullfile(base, 'model_comparisons', ['Pop_Rate_F0_' target '.mat']),...
		% 	"pop_rate_F0")
		load(fullfile(base, 'model_comparisons', ['Model_Pop_Timing_F0_' target '.mat']), ...
			"accur_all","C_all", "num_neurons")
		load(fullfile(base, 'model_comparisons', ['Model_Pop_Rate_F0_' target '.mat']),...
			"pop_rate_F0")
	elseif iinstru == 2
		% load(fullfile(base, 'model_comparisons', 'pop_timing_F0_Oboe_subset2.mat'), ...
		% 	"accur_all","C_all", "num_neurons", "mean_acc", "std_acc")
		% accur_all2 = accur_all;
		% C_all2 = C_all;
		% num_neurons2 = num_neurons;
		% mean_acc2 = mean_acc;
		% std_acc2 = std_acc;
		% load(fullfile(base, 'model_comparisons', ['pop_timing_F0_' target '_subset.mat']), ...
		% 	"accur_all","C_all", "num_neurons")
		% load(fullfile(base, 'model_comparisons', ['Pop_Rate_F0_' target '.mat']),...
		% 	"pop_rate_F0")
		load(fullfile(base, 'model_comparisons', ['Model_Pop_Timing_F0_' target '.mat']), ...
			"accur_all","C_all", "num_neurons")
		load(fullfile(base, 'model_comparisons', ['Model_Pop_Rate_F0_' target '.mat']),...
			"pop_rate_F0")
	elseif iinstru == 3
		% load(fullfile(base, 'model_comparisons', 'pop_timing_F0_invariant.mat'), ...
		% 	"accur_all","C_all", "num_neurons")
		load(fullfile(base, 'model_comparisons', 'Model_Pop_Timing_F0_Invariant.mat'), ...
			"accur_all","C_all", "num_neurons")
	end


	% Get stimulus
	if iinstru == 1 || iinstru == 2
		F0s = getF0s(target);
	else
		F0s = getF0s('Bassoon');
		F0s = F0s(25:40);
	end

	%% Plot

	% if iinstru == 1 || iinstru == 3
		nmodels = length(num_neurons);
		mean_acc = mean(accur_all,2);
		std_acc = std(accur_all, [], 2);
	% else
	% 	nmodels = length(num_neurons) + length(num_neurons2);
	% 	mean_acc1 = mean(accur_all,2);
	% 	mean_acc = [mean_acc1(1:12); mean_acc2(1:3); mean_acc1(13:end); mean_acc2(4:6)];
	% 	std_acc1 = std(accur_all, [], 2);
	% 	std_acc = [std_acc1(1:12); std_acc2(1:3); std_acc1(13:end); std_acc2(4:6)];
	% 	num_neurons = [num_neurons(1:12) num_neurons2(1:3) num_neurons(13:24) num_neurons2(4:6)];
	% end

	h(hind(iinstru)) = subplot(3, 2, hind(iinstru));
	errorbar(num_neurons(1:nmodels/2), mean_acc(1:nmodels/2), std_acc(1:nmodels/2)/sqrt(10));
	hold on
	errorbar(num_neurons(nmodels/2+1:end)', mean_acc(nmodels/2+1:end), std_acc(nmodels/2+1:end)/sqrt(10))
	xlabel('# Neurons in Model')
	ylabel('Accuracy')
	if iinstru == 1
		title('Model Accuracy')
	end
	ylim([0 1])
	grid on
	set(gca, 'FontSize', fontsize)
	box off
	if iinstru == 1
		xlim([0 40])
	else
		xlim([0 100])
	end

	%%

	C_diag_all = NaN(nmodels, length(F0s));
	if iinstru == 1 || iinstru == 2
		for ii = 1:nmodels
			C_diag = NaN(1, length(F0s));
			for irep = 1 %:10
				C_diag(irep,:) = diag(C_all{ii,irep});
			end

			C_diag_all(ii,:) = C_diag; %mean(C_diag);
		end
	% elseif iinstru == 2
	% 	iii = 1;
	% 	iiii = 1;
	% 	for ii = 1:nmodels
	% 		if ~ismember(ii, [13, 14, 15, 28, 29, 30])
	% 			C_diag = NaN(10, length(F0s));
	% 			for irep = 1:10
	% 				C_diag(irep,:) = diag(C_all{iiii,irep});
	% 			end
	% 			iiii = iiii +1;
	% 			C_diag_all(ii,:) = mean(C_diag);
	% 		else
	% 			C_diag_all(ii,:) = diag(C_all2{iii});
	% 			iii = iii +1;
	% 		end
	% 	end
	else
		for ii = 1:nmodels
			C_diag_all(ii,:) = diag(C_all{ii});
		end
	end
	y = F0s;
	x = num_neurons(1:nmodels/2);

	h(hind(iinstru)+1) = subplot(3, 2, hind(iinstru)+1);
	pcolor(y, x, C_diag_all(1:nmodels/2,:), 'EdgeColor','none')
	set(gca, 'xscale', 'log')
	if iinstru == 3
		clim([0 40])
	else
		clim([0 20])
	end
	ylabel('# Neurons in Model')
	%shading interp
	xlabel('F0 (Hz)')
	if iinstru == 1
		title('Best Units')
	end
	a=colorbar;
	a.Label.String = '# Accurate Predictions';
	xticks([50 100 250 500 1000 1500])
	xticklabels([50 100 250 500 1000 1500]/1000)
	set(gca, 'FontSize', fontsize)
	box off

	%%
	% 
	% h(hind(iinstru)+2) = subplot(2, 3, hind(iinstru)+2);
	% pcolor(y, x, C_diag_all(nmodels/2+1:end,:), 'EdgeColor','none')
	% set(gca, 'xscale', 'log')
	% clim([0 20])
	% ylabel('# Neurons in Model')
	% xlabel('F0 (Hz)')
	% title('Worst Units')
	% a=colorbar;
	% a.Label.String = '# Accurate Predictions';
	% xticks([50 100 250 500 1000 1500])
	% xticklabels([50 100 250 500 1000 1500]/1000)
	% set(gca, 'FontSize', fontsize)
	% box off

	%% Extra
	% figure
	% C_timing = C_all{12, 1};
	% response = pop_rate_F0.response;
	% order = F0s;
	% 
	% pred_time = zeros(size(response)); % Preallocate
	% idx = 1;
	% for i = 1:length(F0s)
	% 	true_class = order(i);
	% 	% Find indices of samples with this true class
	% 	true_idx = find(response == true_class);
	% 	% For each predicted class, assign the appropriate number of predictions
	% 	count = 0;
	% 	for j = 1:length(F0s)
	% 		pred_class = order(j);
	% 		n = C_timing(i,j);
	% 		if n > 0
	% 			pred_time(true_idx(count+1:count+n)) = pred_class;
	% 			count = count + n;
	% 		end
	% 	end
	% end
	% pred_rate = pop_rate_F0.validationPredictions;
	% pred_rate2 = response' - pred_rate;
	% pred_time2 = response - pred_time;
	% 
	% scatter(response, pred_time)
	% hold on
	% scatter(response, pred_rate)
	% set(gca, 'xscale', 'log', 'yscale', 'log')


end

%% Arrange positions

left = [0.15 0.58]; 
bottom = linspace(0.08, 0.73, 3);
height = 0.2;
width = 0.25;

set(h(1), 'position', [left(1) bottom(3) width height])
set(h(2), 'position', [left(2) bottom(3) width height])
set(h(3), 'position', [left(1) bottom(2) width height])
set(h(4), 'position', [left(2) bottom(2) width height])
set(h(5), 'position', [left(1) bottom(1) width height])
set(h(6), 'position', [left(2) bottom(1) width height])


%% Annotate

annotation("textbox", [0.11 0.7 0.1386 0.1088], "String", "Bassoon",...
	"FontSize", 12, "FontWeight", "bold", "EdgeColor", "none", "Rotation",90)
annotation("textbox", [0.055 0.43 0.09778 0.05666], "String", "Oboe",...
	"FontSize", 12, "FontWeight", "bold", "EdgeColor", "none", "Rotation",90)
annotation("textbox", [0.11 0.11 0.1386 0.1088], "String", "Both",...
	"FontSize", 12, "FontWeight", "bold", "EdgeColor", "none", "Rotation",90)

labelleft = left-0.13;
annotation('textbox',[labelleft(1) 0.95 0.071 0.058],...
	'String','A','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(2) 0.95 0.071 0.058],...
	'String','B','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
% annotation('textbox',[labelleft(3) 0.96 0.071 0.058],...
% 	'String','C','FontWeight','bold','FontSize',labelsize,...
% 	'EdgeColor','none');


%% Save figure

if save_fig == 1
	filename = 'fig9_pitch_population_timing';
	save_figure(filename)
end
%% supp_image_plots


%% Load in data

% Load in oboe
[base, datapath, ~, ppi] = getPathsNT();

% load(fullfile(base, 'model_comparisons', 'Neuron_Time_F0_Oboe.mat'), ...
% 	"neuron_time_F0")
% neuron_time_F0_oboe = neuron_time_F0;
% load(fullfile(base, 'model_comparisons','Neuron_Rate_F0_Oboe.mat'), ...
% 	"neuron_rate_F0")
% neuron_rate_F0_oboe = neuron_rate_F0;
% 
% % Load in bassoon
% load(fullfile(base, 'model_comparisons', 'Neuron_Time_F0_Bassoon.mat'), ...
% 	"neuron_time_F0")
% load(fullfile(base, 'model_comparisons','Neuron_Rate_F0_Bassoon.mat'), ...
% 	"neuron_rate_F0")
% load(fullfile(base,'model_comparisons', 'Data_NT_3.mat'), 'nat_data')

% Model
load(fullfile(base, 'model_comparisons', 'Model_Neuron_Time_F0_Oboe.mat'), ...
	"neuron_time_F0")
neuron_time_F0_oboe = neuron_time_F0;
load(fullfile(base, 'model_comparisons','Model_N_Rate_F0_Oboe.mat'), ...
	"neuron_rate_F0")
neuron_rate_F0_oboe = neuron_rate_F0;

% Load in bassoon
load(fullfile(base, 'model_comparisons', 'Model_Neuron_Time_F0_Bassoon.mat'), ...
	"neuron_time_F0")
load(fullfile(base, 'model_comparisons','Model_N_Rate_F0_Bassoon.mat'), ...
	"neuron_rate_F0")
load(fullfile(base,'model_comparisons', 'Data_NT_3.mat'), 'nat_data')

%% Set up figure

figure('Position',[50 50 ppi*6, ppi*4])
fontsize = 8;

% Plot
accuracy_rate = [neuron_rate_F0.accuracy]*100;
accuracy_time = [neuron_time_F0.accuracy]*100;
accuracy_rate_oboe = [neuron_rate_F0_oboe.accuracy]*100;
accuracy_time_oboe = [neuron_time_F0_oboe(1:183).accuracy]*100;

for ineuron = 1:4

	if ineuron == 1
		target = 'Bassoon';
		type = 'Timing';
		[~,best_ind] = sort(abs(accuracy_time), 'ascend' );
		C_acc = NaN(length(best_ind), 40);
		for ii = 1:length(best_ind)
			C_acc(ii,:) = diag(neuron_time_F0(best_ind(ii)).C);
		end
	elseif ineuron == 2
		target = 'Bassoon';
		type = 'Rate';
		[~,best_ind] = sort(abs(accuracy_rate), 'ascend' );
		C_acc = NaN(length(best_ind), 40);
		for ii = 1:length(best_ind)
			C_acc(ii,:) = diag(neuron_rate_F0(best_ind(ii)).C);
		end
	elseif ineuron == 3
		target = 'Oboe';
		type = 'Timing';
		[~,best_ind] = sort(abs(accuracy_time_oboe), 'ascend' );
		C_acc = NaN(length(best_ind), 35);
		for ii = 1:length(best_ind)
			C_acc(ii,:) = diag(neuron_time_F0_oboe(best_ind(ii)).C);
		end
	else
		target = 'Oboe';
		type = 'Rate';
		[~,best_ind] = sort(abs(accuracy_rate_oboe), 'ascend' );
		C_acc = NaN(length(best_ind), 35);
		for ii = 1:length(best_ind)
			C_acc(ii,:) = diag(neuron_rate_F0_oboe(best_ind(ii)).C);
		end
	end

	h(ineuron) = subplot(1, 4, ineuron);
	x = getF0s(target);
	y = 1:size(C_acc, 1);

	pcolor(x, y, C_acc, 'EdgeColor','none', 'EdgeAlpha',0)
	set(gca, 'xscale', 'log')
	xticks([60 100 200 350 550 1000 1500])
	xlabel('F0 (Hz)')
	ylabel('Neurons Sorted by Accuracy')
	%c = colorbar;
	title([target ' ' type])
	set(gca, 'fontsize', fontsize)
	%c.Label.String = '# Accurate Predictions';

end

%% Annotate

bottom = 0.1;
left = linspace(0.08, 0.82, 4);
width = 0.17;
height = 0.85;
set(h(1), 'position', [left(1) bottom width height])
set(h(2), 'position', [left(2) bottom width height])
set(h(3), 'position', [left(3) bottom width height])
set(h(4), 'position', [left(4) bottom width height])


labelleft= linspace(0.03, 0.78, 4);
annotation('textbox',[labelleft(1) 0.95 0.071 0.058],...
	'String','A','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(2) 0.95 0.071 0.058],...
	'String','B','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(3) 0.95 0.071 0.058],...
	'String','C','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(4) 0.95 0.071 0.058],...
	'String','D','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');

%% Save figure

save_fig = 0;
if save_fig == 1
	filename = 'supp3_image_acc_plots';
	save_figure(filename)
end
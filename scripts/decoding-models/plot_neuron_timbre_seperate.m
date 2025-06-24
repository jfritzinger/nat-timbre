%% plot_neuron_timbre_seperate
clear
[base, datapath, savepath, ppi] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Neuron_Rate_Timbre_Separate.mat'), ...
	"neuron_rate_timbre")
load(fullfile(base, 'model_comparisons', 'Neuron_Time_Timbre_Separate.mat'), ...
	"neuron_time_timbre")

%% Arrange into matrix

for iii = 1:246
	for ii = 1:16
		acc_matrix(iii, ii) = neuron_rate_timbre(iii, ii).accuracy;
	end
end

figure
tiledlayout(1, 2)
nexttile
pcolor(1:16, 1:246, acc_matrix, 'EdgeColor','none');
colorbar
title('Rate Model')
xlabel('F0s')
ylabel('Neurons')
clim([0 1])

% %% Sort by CF 
% 
% CFs = [neuron_rate_timbre(:,1).CF];
% [CFs, order] = sort(CFs);
% 
% acc_matrix_sorted = acc_matrix(order,:);
% 
% figure
% % pcolor(1:16, CFs, acc_matrix_sorted, 'EdgeColor','none');
% % set(gca, 'YScale', 'log')
% pcolor(1:16, 1:246, acc_matrix_sorted, 'EdgeColor','none');


for iii = 1:246
	for ii = 1:16
		acc_matrix_time(iii, ii) = neuron_time_timbre(iii, ii).accuracy;
	end
end

nexttile
pcolor(1:16, 1:246, acc_matrix_time, 'EdgeColor','none');
colorbar
title('Timing Model')
xlabel('F0s')
ylabel('Neurons')
clim([0 1])

%% 

figure
tiledlayout(4, 4)

for ii = 1:16

	acc_rate = acc_matrix(:,ii);
	acc_time = acc_matrix_time(:,ii);

	nexttile
	scatter(acc_rate, acc_time, 10, 'filled')
	hold on
	plot([0 1], [0 1])
	xlim([0 1])
	ylim([0 1])
	xlabel('Rate accuracy')
	ylabel('Timing accuracy')
end
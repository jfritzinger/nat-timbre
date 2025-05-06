%% plot_neuron_time
clear

% Currently rerunning this on H2 to get accurate models 
%% Load in data 
target = 'Oboe';

[base, datapath, savepath, ppi] = getPathsNT();
load(fullfile(base, 'model_comparisons', ['Neuron_Time_F0_' target '.mat']), ...
	"neuron_time_F0")

%% Plot accuracy of each neuron
figure('Position',[560,618,798,230])
tiledlayout(1, 3)
accuracy(1,:) = [neuron_time_F0.accuracy]*100;
accuracy(2,:) = [neuron_time_F0.accuracy_low]*100;
accuracy(3,:) = [neuron_time_F0.accuracy_high]*100;

for ii = 1:3

	nexttile
	edges = linspace(0, 100, 30);
	histogram(accuracy(ii,:),edges)
	hold on
	xline(2.5, 'k')
	xline(mean(accuracy(ii,:)), 'r', 'LineWidth',2)
	ylabel('# Neurons')
	xlabel('Prediction Accuracy (%)')
	title('Prediction of F0')
	xlim([0 100])

end

%% Find indices of the best neurons 

[temp,originalpos] = sort(validationAccuracy, 'descend' );
n = temp(1:3);
best_ind=originalpos(1:3);

putatives = nat_data(best_ind(1)).putative;
figure
confusionchart(C)
title(num2str(accur(ind)*100))


%% Plot rasters of three best

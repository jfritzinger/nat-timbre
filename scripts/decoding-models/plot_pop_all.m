%% plot_pop_all
clear 

%% Load in data 

base = getPathsNT;
load(fullfile(base, 'model_comparisons', 'Pop_Rate.mat'), "pop_rate")


%% Confusion matrix

% Compute classification using all data 
figure
confusionchart(pop_rate.C)

% Calculate accuracy
chart = confusionchart(pop_rate.response,pop_rate.validationPredictions); % Generate confusion chart
accuracy(1) = sum(diag(pop_rate.C)) / sum(pop_rate.C(:)); % Calculate accuracy
%title(sprintf('Accuracy = %0.2f%%', accuracy(1)*100))

%% Load in timing data 

base = getPathsNT;
load(fullfile(base, 'model_comparisons', 'Pop_Time.mat'), ...
	"accur_all","C_all", "num_neurons", "mean_acc", "std_acc")

figure
tiledlayout(1, 2)
nexttile 
plot(num_neurons, accur_all)
ylabel('Overall Accuracy')
xlabel('# Neurons in Model')
title('Model using populations of PSTHs')

nexttile

for ind = 1:10
	C_diag(ind,:) = diag(C_all{ind});
end
pcolor(num_neurons, 1:75,C_diag', 'EdgeColor','none')
hold on
yline(40, 'LineWidth',3)
colorbar
title('Correct predictions for each category')
ylabel('Category (Instrument & F0)')
xlabel('# Neurons in model')

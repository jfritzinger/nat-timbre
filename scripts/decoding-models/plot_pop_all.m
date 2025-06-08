%% plot_pop_all

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
title(sprintf('Accuracy = %0.2f%%', accuracy(1)*100))
%% plot_pop_rate_F0
clear

%% Load data 
target = 'Oboe';

[base, ~, ~, ~] = getPathsNT();
load(fullfile(base, 'model_comparisons', ['Pop_Rate_F0_' target '.mat']),...
	"pop_rate_F0")

%% Plot many model reps 
% 
% figure
% swarmchart(ones(50,1), validationAccuracy, 'filled')
% hold on
% boxplot(validationAccuracy)
% ylabel('Model Accuracy')
% title('50 reps of model')

%% Plot confusion matrix 

% Compute classification using all data 
figure
confusionchart(pop_rate_F0.C)

% Calculate accuracy
chart = confusionchart(pop_rate_F0.response,pop_rate_F0.validationPredictions); % Generate confusion chart
confusionMatrix = chart.NormalizedValues; % Get the normalized confusion matrix
accuracy(1) = sum(diag(confusionMatrix)) / sum(confusionMatrix(:)); % Calculate accuracy
title(sprintf('Accuracy = %0.2f%%', accuracy(1)*100))

% Accurate % 202 / 800
num_acc =  sum(diag(confusionMatrix));

% Accuracy within +/- 1 semitone % 13 / 800
num_acc1 =  sum(diag(confusionMatrix, 1)) + sum(diag(confusionMatrix, -1));

% Accuracy within +/- 2 semitones % 319 / 800
num_acc2 =  sum(diag(confusionMatrix, 2)) + sum(diag(confusionMatrix, -2));

% Accuracy within +/- 3 semitones % 319 / 800
num_acc3 =  sum(diag(confusionMatrix, 3)) + sum(diag(confusionMatrix, -3));

% Accuracy within +/- 4 semitones % 319 / 800
num_acc4 =  sum(diag(confusionMatrix, 4)) + sum(diag(confusionMatrix, -4));


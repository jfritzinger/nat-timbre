%% 
clear

%% Load in data

filepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/model_comparisons';
load(fullfile(filepath, 'Data_NT.mat'), 'nat_data')

%% Get correct output of model

if ismac
	fpath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data';
else
	fpath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\';
end
tuning = readtable(fullfile(fpath, 'Tuning.xlsx')); % Load in tuning

% Get bassoon stimulus
target = 'Bassoon';
listing = dir(fullfile(fpath, 'waveforms', ['*' target '*.wav']));
files = {listing.name};
note_names = extractBetween(files, 'ff.', '.');
[~, index] = ismember(note_names, tuning.Note);
F0s1 = tuning.Frequency(index);
[F0s1, order] = sort(F0s1);
response = reshape(repmat(F0s1, 1, 20)', 1, []);
response = log10(response)';

%% Get all rates for each repetition for bassoon (one example neuron)


putative = 'R29_TT3_P2_N06'; % Good
index = find(strcmp({nat_data.putative}, putative));

% Get data
h_all = [];
for itarget = 1:40
	spikes_bass = nat_data(index).bass_spikerate{itarget}/1000; % ms
	spikereps_bass = nat_data(index).bass_spikerep{itarget};

	% Arrange data for SVM
	min_dis = 0.25;
	edges = 0:min_dis:300;
	t = 0+min_dis/2:min_dis:300-min_dis/2;
	for irep = 1:20
		h_bass(irep, :) = histcounts(spikes_bass(spikereps_bass==irep), edges);
	end
	h_all = [h_all; h_bass];
end

% Put data into table
T = array2table(h_all);
T.response = response;


%% Split into training and testing data 

for group = 0:39

    % Calculate group range (1-20, 21-40, etc.)
    startIdx = group*20 + 1;
    endIdx = (group+1)*20;
    
    % Random permutation within current group
    shuffled = randperm(20) + startIdx - 1;
    
    % Store indices
    trainStart = group*16 + 1;
    trainEnd = (group+1)*16;
    testStart = group*4 + 1;
    testEnd = (group+1)*4;
    
    trainIndices(trainStart:trainEnd) = shuffled(1:16);
    testIndices(testStart:testEnd) = shuffled(17:20);
end

% Create training and testing sets
T_train = T(trainIndices, :);
T_test = T(testIndices, :);

%% Fit using fitrlinear 

% hyperopts = struct('AcquisitionFunctionName','expected-improvement-plus', 'KFold',5);
% [Mdl,FitInfo,HyperparameterOptimizationResults] = fitrlinear(T_test, 'response', ...
%     'OptimizeHyperparameters','auto', ...
%     'HyperparameterOptimizationOptions',hyperopts);

Mdl = fitrlinear(T, 'response','BetaTolerance',0.0001, ...
	'Learner','leastsquares', 'Lambda','auto', 'Solver','lbfgs', ...
	'KFold',5, 'CrossVal','on', 'Regularization','ridge');

%% Save

pred_F0 = kfoldPredict(Mdl);
%pred_F0 = predict(Mdl, T_test);
[base, datapath, ~, ppi] = getPathsNT();
save(fullfile(base, 'model_comparisons', 'Best_Linear_F0.mat'), "pred_F0", "T")

%% 

r = corrcoef(pred_F0, T.response);
r2 = r(1, 2)^2;

figure('Position',[31,910,910,413])
nexttile
scatter(T.response, pred_F0, 'filled', 'MarkerFaceAlpha',0.5)
%set(gca, 'xscale', 'log', 'yscale', 'log')
% ylim([57 588])
% xlim([57 588])
hold on
plot(T_test.response, T_test.response)
title(['R^2 = ' num2str(r2)])
ylabel('Predicted log10(F0)')
xlabel('Actual log10(F0)')

logF0 = log10(F0s1);
for ii = 1:length(pred_F0)

	differences = abs(logF0 - pred_F0(ii));
	[~, closest_column_index] = min(differences);
	closest_cat(ii) = logF0(closest_column_index);
end


% Plot confusion matrix
nexttile
C = confusionmat(T.response, closest_cat);
confusionchart(C)

% Calculate accuracy
chart = confusionchart(T.response,closest_cat); % Generate confusion chart
confusionMatrix = chart.NormalizedValues; % Get the normalized confusion matrix
accuracy = sum(diag(confusionMatrix)) / sum(confusionMatrix(:)); % Calculate accuracy
title(sprintf('Accuracy = %0.2f%%', accuracy*100))

